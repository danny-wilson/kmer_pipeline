#!/usr/bin/env python3
"""Translate kmer_pipeline's Dockerfile into an Apptainer definition file.

Lets the release image be pre-tested with Apptainer, on a host that has no
Docker. Handles only the instructions the Dockerfile uses: FROM (by digest),
ARG, LABEL, USER, WORKDIR, RUN, COPY . <dir>, ENV. RUN commands run in order
with /bin/sh -c from the current WORKDIR, as Docker does. COPY is done in
%files, before all of %post, so the script checks that no earlier RUN refers
to the COPY destination.

Usage: dockerfile2def.py DOCKERFILE CONTEXT_DIR VERSION > image.def
"""
import re
import shlex
import sys


def logical_lines(text):
    out, cur = [], ""
    for line in text.splitlines():
        if not cur and (not line.strip() or line.lstrip().startswith("#")):
            continue
        if line.endswith("\\"):
            cur += line[:-1]
        else:
            out.append(cur + line)
            cur = ""
    assert not cur, "dangling continuation"
    return out


def main(dockerfile, context, version):
    args = {"VERSION": version}
    header, labels, files, post, env = [], [], [], [], []
    workdir = "/"
    for line in logical_lines(open(dockerfile).read()):
        instr, _, rest = line.partition(" ")
        rest = rest.strip()
        sub = lambda s: re.sub(r"\$\{(\w+)\}", lambda m: args.get(m.group(1), m.group(0)), s)
        if instr == "FROM":
            header = ["Bootstrap: docker", "From: " + rest]
        elif instr == "ARG":
            assert rest in args, rest
        elif instr == "LABEL":
            k, v = rest.split("=", 1)
            labels.append(f"    {k} {sub(v).strip(chr(34))}")
        elif instr == "USER":
            pass  # Apptainer runs as the caller
        elif instr == "WORKDIR":
            workdir = rest
        elif instr == "RUN":
            assert "'" not in rest, "RUN with a single quote"
            post.append(f"    mkdir -p {workdir} && cd {workdir} && /bin/sh -c '{rest}' || exit 1")
        elif instr == "COPY":
            src, dest = rest.split()
            assert src == ".", rest
            # %files runs before %post: no RUN may come before this COPY and use dest
            assert not any(dest in p for p in post), "COPY destination used by an earlier RUN"
            files.append(f"    {context} {dest}")
        elif instr == "ENV":
            if "=" in rest.split()[0]:
                k, v = rest.split("=", 1)
            else:
                k, v = rest.split(None, 1)
            if k == "HOME":
                continue  # Apptainer sets HOME itself
            env.append(f"    export {k}={shlex.quote(v)}")
            # Docker's ENV also applies to later RUN steps
            post.append(f"    export {k}={shlex.quote(v)}")
        else:
            sys.exit(f"unsupported instruction: {instr}")
    print("\n".join(header))
    print("\n%labels\n" + "\n".join(labels))
    print("\n%files\n" + "\n".join(files))
    print("\n%post\n" + "\n".join(post))
    print("\n%environment\n" + "\n".join(env))


if __name__ == "__main__":
    main(*sys.argv[1:])
