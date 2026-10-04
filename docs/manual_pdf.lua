-- pandoc filter for building the PDF manual from docs/manual.md (read as gfm).
-- The manual uses a little HTML that GitHub renders but LaTeX output would drop:
-- <sub>, <sup>, <em>, <br>, and centred images and captions written as
--   <p align="center"><img src="..." width="50%" alt="..."></p>
--   <p align="center"><em>Caption.</em></p>
-- This filter turns them into pandoc elements; any other HTML stops the build.
-- It also prints each one-cell "Output" table as an "Output" block between rules (a table row
-- cannot break across pages), and marks every code block as bash so that LaTeX wraps its
-- long lines (manual_pdf.tex). Sections in the Contents start on new pages, and horizontal
-- rules span the full text width.

local function fail(msg)
  io.stderr:write("manual_pdf.lua: " .. msg .. "\n")
  os.exit(1)
end

local spans = { sub = pandoc.Subscript, sup = pandoc.Superscript, em = pandoc.Emph }

-- Horizontal rules span the full text width (pandoc's default is half width)
local function full_rule()
  return pandoc.RawBlock("latex", "\\par\\noindent\\rule{\\linewidth}{0.4pt}\\par")
end

-- Inline HTML: <br> and paired <sub>/<sup>/<em> tags
function Inlines(inlines)
  local out, stack = pandoc.List(), {}
  for _, el in ipairs(inlines) do
    local tag = el.t == "RawInline" and el.format == "html" and el.text
    if tag == "<br>" then
      (#stack > 0 and stack[#stack].content or out):insert(pandoc.LineBreak())
    elseif tag and tag:match("^<(%a+)>$") and spans[tag:match("^<(%a+)>$")] then
      table.insert(stack, { name = tag:match("^<(%a+)>$"), content = pandoc.List() })
    elseif tag and tag:match("^</(%a+)>$") then
      local top = table.remove(stack)
      if not top or top.name ~= tag:match("^</(%a+)>$") then fail("unbalanced " .. tag) end
      local made = spans[top.name](top.content);
      (#stack > 0 and stack[#stack].content or out):insert(made)
    elseif tag then
      fail("HTML not handled: " .. tag)
    else
      (#stack > 0 and stack[#stack].content or out):insert(el)
    end
  end
  if #stack > 0 then fail("unclosed <" .. stack[#stack].name .. ">") end
  return out
end

-- Block HTML: a centred image followed by its centred caption becomes a figure
function Blocks(blocks)
  local out, i = pandoc.List(), 1
  while i <= #blocks do
    local b = blocks[i]
    if b.t == "RawBlock" and b.format == "html" then
      local src, width, alt = b.text:match('^<p align="center"><img src="([^"]+)" width="([^"]+)" alt="([^"]*)"></p>%s*$')
      if src then
        local nxt = blocks[i + 1]
        local cap = nxt and nxt.t == "RawBlock" and nxt.format == "html"
          and nxt.text:match('^<p align="center"><em>(.-)</em></p>%s*$')
        if not cap then fail("image without a centred caption: " .. src) end
        local img = pandoc.Image(pandoc.read(cap, "gfm").blocks[1].content, src, "fig:" .. alt)
        img.attributes.width = "75%" -- half width on GitHub is too small on a page
        out:insert(pandoc.Para({ img }))
        i = i + 2
      else
        fail("HTML block not handled: " .. b.text:sub(1, 60))
      end
    else
      out:insert(b)
      i = i + 1
    end
  end
  return out
end

function CodeBlock(el)
  el.classes = { "bash" }
  return el
end

-- Column widths in proportion to each column's longest cell (capped), so LaTeX wraps the
-- cells instead of running tables off the page (gfm tables have no widths)
local function set_widths(el)
  local rows = pandoc.List()
  rows:extend(el.head.rows)
  for _, b in ipairs(el.bodies) do rows:extend(b.body) end
  local len = {}
  for _, row in ipairs(rows) do
    for i, cell in ipairs(row.cells) do
      local n = math.min(math.max(#pandoc.utils.stringify(cell.contents), 6), 70)
      len[i] = math.max(len[i] or 6, n)
    end
  end
  local total = 0
  for _, n in ipairs(len) do total = total + n end
  for i, spec in ipairs(el.colspecs) do
    el.colspecs[i] = { spec[1], 0.97 * len[i] / total }
  end
  return el
end

-- One-cell "Output" table: rule, bold "Output", the cell's text split into paragraphs, rule
function Table(el)
  local head = el.head.rows
  local body = el.bodies[1] and el.bodies[1].body or {}
  if #head == 1 and #head[1].cells == 1 and #body == 1 and #body[1].cells == 1
      and pandoc.utils.stringify(head[1].cells[1].contents) == "Output" then
    local out = pandoc.List({ full_rule(), pandoc.Para({ pandoc.Strong("Output") }) })
    for _, b in ipairs(body[1].cells[1].contents) do
      local para = pandoc.List()
      local breaks = 0
      for _, x in ipairs(b.content) do
        if x.t == "LineBreak" then
          breaks = breaks + 1
        else
          if breaks > 0 and #para > 0 then
            out:insert(pandoc.Para(para))
            para = pandoc.List()
          end
          breaks = 0
          para:insert(x)
        end
      end
      if #para > 0 then out:insert(pandoc.Para(para)) end
    end
    out:insert(full_rule())
    return out
  end
  if #el.colspecs > 6 then
    -- wide tables (the k-mer table): natural column widths, scaled to the text width
    local function tex(x)
      return (pandoc.utils.stringify(x):gsub("[\\{}%%$&#_^~]", function(c) return "\\" .. c .. "{}" end))
    end
    local lines = {}
    local function row(cells)
      local t = {}
      for _, c in ipairs(cells) do table.insert(t, tex(c.contents)) end
      table.insert(lines, table.concat(t, " & ") .. " \\\\")
    end
    for _, r in ipairs(el.head.rows) do row(r.cells) end
    table.insert(lines, "\\hline")
    for _, b in ipairs(el.bodies) do for _, r in ipairs(b.body) do row(r.cells) end end
    local spec = "@{}l" .. string.rep("r", #el.colspecs - 1) .. "@{}"
    return pandoc.RawBlock("latex", "\\begin{center}\\resizebox{\\textwidth}{!}{\\begin{tabular}{" .. spec
      .. "}\\hline\n" .. table.concat(lines, "\n") .. "\n\\hline\\end{tabular}}\\end{center}")
  end
  return set_widths(el)
end

-- Each section listed in the Contents (level-2 headings) starts on a new page
function Header(el)
  if el.level == 2 and pandoc.utils.stringify(el) ~= "Contents" then
    return { pandoc.RawBlock("latex", "\\clearpage"), el }
  end
end

function HorizontalRule()
  return full_rule()
end
