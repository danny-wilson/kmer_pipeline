// Based on https://www.w3schools.com/howto/howto_js_slideshow.asp
// and https://stackoverflow.com/questions/43299759/how-do-i-make-multiple-slideshows-in-the-same-html-document
// Get all slideshow containers and record a slideIndex for each
var slideShows = document.getElementsByClassName("slideshow-container")
//console.log("Number of slideshows:")
//console.log(slideShows.length)
//let slideIndex = new Array()
let j
for(j=0;j<slideShows.length;j++) {
	slideShows[j].slideIndex = 1
}
for(j=0;j<slideShows.length;j++) {
	showSlides(slideShows[j].slideIndex, slideShows[j])
}
showToggled()

//let slideIndex = 1;
//showSlides(slideIndex);
function plusSlides(n,elem) {
	//console.log(elem)
	slideShow = elem.closest(".slideshow-container")
	//console.log("plusSlides(): elem = " + elem.class)
	showSlides(slideShow.slideIndex += n, slideShow);
}
function currentSlide(n,elem) {
  console.log(elem)
  slideShow = elem.closest(".slideshow-container")
  showSlides(slideShow.slideIndex = n, slideShow);
}
function showSlides(n, slideShow) {
  let i;
  let slides = slideShow.getElementsByClassName("mySlides");
  let dots = slideShow.getElementsByClassName("dot");
  if (n > slides.length) {slideShow.slideIndex = 1}
  if (n < 1) {slideShow.slideIndex = slides.length}
  for (i = 0; i < slides.length; i++) {
    slides[i].style.display = "none";
  }
  for (i = 0; i < dots.length; i++) {
    dots[i].className = dots[i].className.replace(" active", "");
  }
  slides[slideShow.slideIndex-1].style.display = "block";
  dots[slideShow.slideIndex-1].className += " active";
//  showToggled()
}
function showToggled() {
	let toggler = document.getElementById("toggler");
	let checked = document.getElementsByClassName("toggled");
	let unchecked = document.getElementsByClassName("untoggled");
	if(toggler.checked) {
		for(i=0;i<checked.length;i++) checked[i].style.display = "block";
		for(i=0;i<unchecked.length;i++) unchecked[i].style.display = "none";
	} else {
		for(i=0;i<checked.length;i++) unchecked[i].style.display = "block";
		for(i=0;i<unchecked.length;i++) checked[i].style.display = "none";
	}
}
