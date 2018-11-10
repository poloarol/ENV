function plotGenome(data){
  let mainCanvas = document.getElementById("canvas");
  let secondary = document.getElementById("itemMap");
  appendElem(secondary, data);

  document.getElementById("name").innerHTML = data[0][0]["orgName"];

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

  const span = document.querySelectorAll(".ui.compact.icon.button");
  for(let i = 0; i < span.length; i++){
    span[i].addEventListener("click", function(){
      // copyToClickBoard(document.querySelectorAll("ui.compact.icon.button")[i]);
    })
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


function appendElem(container, data){
  for(let i = 1; i < data.length; i++){
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
    let icon = document.createElement("div");
    let geneIcon = document.createElement("i");

    let buttonDiv = document.createElement("button");
    buttonDiv.setAttribute("class", "ui compact icon basic button");
    geneIcon.setAttribute("class", "large cloud download icon");
    // let fastaIcon = document.createElement("i");
    // fastaIcon.setAttribute("class", "cloud download icon");
    buttonDiv.appendChild(geneIcon);
    let geneDiv = document.createElement("div");
    // let fastaDiv = document.createElement("div");
    // let fastaDivText = document.createTextNode("Fasta Seq. ")
    let geneDivText = document.createTextNode("AA Sequence");                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    geneDiv.appendChild(geneDivText);
    geneDiv.appendChild(buttonDiv);
    // fastaDiv.appendChild(fastaDivText);
    // fastaDiv.appendChild(fastaIcon);
    icon.appendChild(geneDiv);
    // icon.appendChild(fastaDiv);
    let name = document.createTextNode(`${data[i][0]["orgName"]}`);
    let accession = document.createTextNode(`Accession #:`);

    let value = sequence(data[i]);
    buttonDiv.setAttribute("data-value", value);

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
    div.appendChild(orgDesc);
    div.appendChild(canvas);
    orgDescContent.appendChild(icon)
    container.appendChild(div);
  }
}

function sequence(data){
  let value = "";
  for(let i = 0; i < data.length; i++){
    value = value + data[i]["AA"]
  }
  return value;
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

/*function copyToClickBoard(element){
  console.log(element);
}*/