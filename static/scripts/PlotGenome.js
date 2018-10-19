function plotGenome(data){
  let mainCanvas = document.getElementById("canvas");

  let array = new Array();
  for(let i = 0; i < 11; i++){
    array.push(data);
  }

  // let secondary = document.getElementById("secondary");
  // appendElem(secondary);

  let chart = new Scribl(mainCanvas, 500);

  for(let i = 0; i < array.length; i++){
    gap_start = find_start_stop(array[i]);
    if(data[gap_start-1] !== null){
      normalize(array[i], gap_start);
    }
  }

  generateDiagram(chart, array.splice(0,1));
  for(let i = 0; i < array.length; i++){
    let canvasID = "canvas" + i;
    let secondaryCanvas = document.getElementById(canvasID);
    let secondaryChart = new Scribl(secondaryCanvas, 500);
    let data = [array[i]]
    generateDiagram(secondaryChart, data);
  }

}

function normalize(data, start){
  let begin = 500;

  for(let i = 0; i < data.length; i++){
    if (i === 0){
      data[i].start = begin;
    }else{
      begin = begin + data[i].length;
      data[i].start = begin;
    }
  }
}


function find_start_stop(data){
  for(let i = 0; i < data.length; i++){
    if(data[i].id === 1){
      return i;
    }
  }
}


// function appendElem(container){
//   for(let i = 1; i < 11; i++){
//     let div = document.createElement("div");
//     let canvas = document.createElement("canvas");
//     let height = document.createAttribute("height");
//     let width = document.createAttribute("width");
//     let canvasID = document.createAttribute("id");
//     height.value = 430;
//     width.value = 650;
//     canvasID.value = "canvas" + i;
//     canvas.setAttributeNode(height);
//     canvas.setAttributeNode(width);
//     canvas.setAttributeNode(canvasID);
//     div.setAttribute("class", "slides");
//     div.appendChild(canvas);
//     container.appendChild(div)
//   }
// }

function generateDiagram(chart, data){
  for(let i = 0; i < 1; i++){
    let track = chart.addTrack();
    for(let j = 0; j < data[i].length; j++){
       let color = (data[i][j].strand === "+") ? "#998ec3" : "#f1a340";
       gene = track.addFeature(new BlockArrow('complex', data[i][j].start, data[i][j].length, data[i][j].strand, {'color': color}));
       gene.onMouseover = "Name: " + data[i][j].name + "";
    }
  }
  chart.draw();
}

let slideIndex = 1;
showSlides(slideIndex);

function plusSlides(n) {
  showSlides(slideIndex += n);
}

function currentSlide(n) {
  showSlides(slideIndex = n);
}

function showSlides(n) {
  let slides = document.getElementsByClassName("slides");
  let dots = document.getElementsByClassName("dot");
  if (n > slides.length) {slideIndex = 1}    
  if (n < 1) {slideIndex = slides.length}
  for (let i = 0; i < slides.length; i++) {
      slides[i].style.display = "none";  
  }
  for (let i = 0; i < dots.length; i++) {
      dots[i].className = dots[i].className.replace(" active", "");
  }
  slides[slideIndex-1].style.display = "block";  
  dots[slideIndex-1].className += " active";
}