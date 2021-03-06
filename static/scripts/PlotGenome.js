function plotGenome(data){
  let mainCanvas = document.getElementById("canvas");
  appendElem(data)
  let slides = createSlides(data);
  for(let i = 0; i < slides.length; i++){
     let slide = "#slides"+i;
     let div = document.querySelectorAll(slide);
     if(div[0].offsetParent.classList.contains("cloned"))
        div[1].appendChild(slides[i]);
     else
        div[0].appendChild(slides[i]);
   }
  document.getElementById("name").innerHTML = data[0][0]["orgName"];
  document.getElementById("acc_num").innerHTML = `Accession number : ${data[0][0]["accession_number"]}`;

  let chart = new Scribl(mainCanvas, 500);

  for(let i = 0; i < data.length; i++){
    gap_start = find_start_stop(data[i]);
    if(data[gap_start-1] !== null){
      normalize(data[i], gap_start);
    }
  }

  generateDiagram(chart, [data[0]]);
  for(let i = 0; i < data.length; i++){
    let canvasID = "#canvas" + i;
    let secondaryCanvas = document.querySelectorAll(canvasID);
    let secondaryChart = new Scribl(secondaryCanvas[0], 500);
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


function createSlides(data){
  // Attach eash slide to at it's position
  let slides = [];
  for(let i = 0; i < data.length; i++){
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
    geneIcon.setAttribute("class", "large copy icon");
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
    let name = document.createTextNode(`${data[i][1]["orgName"]}`);
    let accession = document.createTextNode(`Accession number : ${data[i][0]["accession_number"]}`);

    height.value = 400;
    width.value = 700;

    canvasID.value = "canvas" + i;
    canvas.setAttributeNode(height);
    canvas.setAttributeNode(width);
    canvas.setAttributeNode(canvasID);
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
    div.setAttribute("class", "canvasDiv")
    orgDesc.appendChild(orgDescContent);
    div.appendChild(orgDesc);
    div.appendChild(canvas);
    orgDescContent.appendChild(icon)
    slides.push(div);
  }
  return slides
}

function appendElem(slides){
  // creates empty slides
  let slideID;
  for(let i = 0; i < slides.length; i++){
    slideID = "slides"+i;
    $('#map').owlCarousel()
              .trigger('add.owl.carousel', [jQuery(`<div class="owl-item" id = ${slideID}>` +  + `</div>`)])
              .trigger('refresh.owl.carousel');
  }
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