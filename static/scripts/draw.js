function draw(parsedData){
  var svgContainer = d3.select(".svg").append("svg")
              .attr("height", "100%")
              .attr("width", "100%");

  console.log(parsedData)

  parsedData.sort(function(a, b){
    return a.Start - b.Start
  });


  // Define the x-axis horizontal
  var minValue = parsedData[0].Start;
  var maxValue = 1000 + parsedData[parsedData.length - 1].Stop;

  var xScale = d3.scaleLinear()
                       .domain([minValue, maxValue])
                       .range([50, 1000])

  var line = svgContainer.selectAll("line")
                .data(parsedData)
                .enter().append("svg:line")
                .attr("x1", function(d, i){
                  if(i== 0){
                    return xScale(d.Start)
                  }
                  return xScale(100 + d.Start)})
                .attr("y1", "100")
                .attr("x2", function(d){ return xScale(d.Stop) })
                .attr("y2", "100")
                .attr("stroke", function(d){
                  if(d.Strand === 1){
                    return "#291B2C";
                  }else{
                    return "#CCA969";
                  };
                })
                .attr("stroke-width", "15")
                .attr("class","tooltip")
                .style("pointer-event","all")
                .append("svg:title")
                .text(function(d){
                  return "Gene: " + d.Gene + "\nLocus Tag: " + d.LocusTag + "\nProduct: " + d.Product;
                })


  // tableCreate(parsedData);

  }

  function tableCreate(parsedData){

    var myArray = ["Gene", "Locus tag", "Product"]

    var tableDiv = document.getElementById('table');
    var tbl = document.createElement('table');
    tbl.style.width = "100%";
    tbl.style.border = "1px solid black";
    // var tRow = tbl.insertRow();
    //
    // for (var i = 0; i < myArray.length; i++) {
    //   var tHead = tRow.insertCell();
    //   tHead.appendChild(document.createTextNode(myArray[i]));
    // }

    for (var i = 0; i < parsedData.length; i++) {
      var tr = tbl.insertRow();
      for (var i = 0; i < myArray.length; i++) {
        var th = tr.insertCell();
        if(i == 0){
          th.appendChild(document.createTextNode(parsedData[i].Gene));
        }else if(i == 1){
          th.appendChild(document.createTextNode(parsedData[i].LocusTag));
        }else{
          th.appendChild(document.createTextNode(parsedData[i].Product))
        }
      }
    }

  tableDiv.appendChild(tbl);
  }
