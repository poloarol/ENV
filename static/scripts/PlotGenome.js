function plotGenome(data){
  let mainCanvas = document.getElementById("canvas");

  let secondary = document.getElementById("secondary");
  appendElem(secondary, data.length-1);

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

  genCheckBox(data);

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
    let orgDesc = document.createElement("div");
    let orgDescContent = document.createElement("div");
    let header = document.createElement("div");
    let meta = document.createElement("div");
    let desc = document.createElement("div");

    let name = document.createTextNode("Organism Name: ");
    let accession = document.createTextNode("Accession Number: ");

    height.value = 400;
    width.value = 700;

    canvasID.value = "canvas" + i;
    canvas.setAttributeNode(height);
    canvas.setAttributeNode(width);
    canvas.setAttributeNode(canvasID);
    div.setAttribute("class", "item");
    desc.setAttribute("class", "description");
    desc.setAttribute("class", "sec");
    desc.innerHTML = "Click to see Gene info.";
    meta.setAttribute("class", "meta");
    meta.appendChild(accession);
    header.setAttribute("class", "header");
    header.appendChild(name);
    orgDesc.setAttribute("class", "ui card");
    orgDescContent.setAttribute("class", "content");
    orgDescContent.appendChild(header);
    orgDescContent.appendChild(meta);
    orgDescContent.appendChild(desc);

    orgDesc.appendChild(orgDescContent);

    div.appendChild(canvas);
    div.appendChild(orgDesc);
    container.appendChild(div);
    
  }
}

function genCheckBox(data){
  let div = document.getElementsByClassName("grouped field")[0];
  for(let i = 0; i < data.length; i++){
    let box = document.createElement("div");
    box.setAttribute("class", "ui checkbox");
    let checkbox = document.createElement('input');
    // checkbox.setAttribute("class", "hidden");
    checkbox.setAttribute("tabindex", 0);
    // checkbox.setAttribute("name", "organism");
    let label = document.createElement("label");
    let name = document.createTextNode(`${data[i][0]["orgName"]}`);
    let br = document.createElement("br");
    let value = "";
    checkbox.type = "checkbox"; 
    for(let j = 0; j < data[i].length; j++){
      value += data[i][j]["AA"];
    }
    checkbox.name =  "organism";
    checkbox.value = value;
    label.appendChild(name);
    // <input type="checkbox" tabindex="0" class="hidden">
    // <label>Checkbox</label>
    box.appendChild(checkbox);
    box.appendChild(label);
    box.appendChild(br);
    div.appendChild(box);
  }
}

function getSeq(data){
  // let seq = ""
  // for (key in data){
  //   if (key.hasOwnProperty("AA")){
  //     seq += data["AA"]
  //   }
  // }
  console.log(data);
}

function generateDiagram(chart, data){
  let gene;
  for(let i = 0; i < 1; i++){
    let track = chart.addTrack();
    let name;
    for(let j = 0; j < data[i].length; j++){
       let color = (data[i][j].strand === "+") ? "#998ec3" : "#f1a340";
       gene = track.addFeature(new BlockArrow('complex', data[i][j].start, data[i][j].length, data[i][j].strand, {'color': color}));
       name = data[i][j].name != "N/A" ? data[i][j].name : (data[i][j].extraclass[1] != "N/A" ? data[i][j].extraclass[1] : data[i][j].extraclass[0]);
       gene.onMouseover = name;
       gene.onClick = { 
                          name: data[i][j].name,
                          start: data[i][j].start,
                          length: data[i][j].length,
                          strand: data[i][j].strand,
                          locus: data[i][j].extraclass[1],
                          protein: data[i][j].extraclass[0]
                      };
    }
  }
  chart.draw();
}