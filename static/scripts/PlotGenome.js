function plotGenome(data){
  let mainCanvas = document.getElementById("canvas");

  let secondary = document.getElementById("secondary");
  appendElem(secondary, data.length);

  let chart = new Scribl(mainCanvas, 500);

  for(let i = 0; i < data.length; i++){
    gap_start = find_start_stop(data[i]);
    if(data[gap_start-1] !== null){
      normalize(data[i], gap_start);
    }
  }

  generateDiagram(chart, data.splice(0,1));
  for(let i = 1; i < data.length; i++){
    let canvasID = "canvas" + i;
    let secondaryCanvas = document.getElementById(canvasID);
    let secondaryChart = new Scribl(secondaryCanvas, 500);
    let datum = [data[i]]
    generateDiagram(secondaryChart, datum);
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


function appendElem(container, length){
  for(let i = 1; i < length; i++){
    let div = document.createElement("div");
    let canvas = document.createElement("canvas");
    let height = document.createAttribute("height");
    let width = document.createAttribute("width");
    let canvasID = document.createAttribute("id");
    height.value = 430;
    width.value = 650;
    canvasID.value = "canvas" + i;
    canvas.setAttributeNode(height);
    canvas.setAttributeNode(width);
    canvas.setAttributeNode(canvasID);
    div.setAttribute("class", "item");
    div.appendChild(canvas);
    container.appendChild(div)
  }
}

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