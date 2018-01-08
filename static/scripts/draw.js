function draw(parsedData){
  var svg = d3.select(".svg").append("svg")
              .attr("height", "100%")
              .attr("width", "100%");

  parsedData.sort(function(a, b){
    return a.Start - b.Start
  });

function myPathForward(start, stop){
    [{ "x" : xScale(start + 100), "y": xScale(start + 150)}, // straight left middle
     {"x" : xScale(start + 100), "y": xScale(start + 100)}, // staight left top
     {"x": xScale(stop + 300), "y": xScale(stop + 100)}, // top right
     {"x": xScale(stop + 350), "y": xScale(stop + 150)}, // right pointy point
     {"x": xScale(stop + 300), "y": xScale(stop + 200)}, // bottom right
     {"x": xScale(start + 100), "y": xScale(start + 200)}, // straight bottom left
     {"x": xScale(start + 100), "y": xScale(start + 150)}]; // Back to left
    }

function myPathBackward(start, stop){
    [{"x" : xScale(start + 50), "y": xScale(start + 150) }, // left pointy bit
     {"x" : xScale(start + 100), "y": xScale(start + 100) }, // left top
     {"x": xScale(stop + 300), "y": xScale(stop + 100) }, // straight top right
     {"x": xScale(stop + 300), "y": xScale(stop + 150) }, // right straight point
     {"x": xScale(stop + 300), "y": xScale(stop + 200) }, // bottom right
     {"x": xScale(start + 100), "y": xScale(start + 200) }, // straight bottom left
     {"x": xScale(start+ 50),  "y": xScale(start + 150)}, ]  // Back to left
    }

  // Define the x-axis horizontal
  var minValue = parsedData[0].Start;
  var maxValue = 1000 + parsedData[parsedData.length - 1].Stop;

  var xScale = d3.scaleLinear()
                       .domain([minValue, maxValue])
                       .range([50, 1500])

  var lineFunction = d3.line()
                       .x(function(d) {return d.x })
                       .y(function(d) {return d.y })
                       .curve(d3.curveBasics);

for (var i = 0; i < parsedData.length; i++) {
  console.log(parsedData[i]);
}

  var line = svg.selectAll("path")
                .data(parsedData)
                .enter().append("svg:path")
                .attr("stroke", function(d){
                  if(d.Strand == 1){
                    return "#291B2C";
                  }else{
                    return "#CCA969";
                  };
                })
                .attr("stroke-width", 1.2)
                .attr("opacity", 0.5)
                .attr("fill", function(d){
                  if(d.Strand == 1){
                    return "#291B2C";
                  }else{
                    return "#CCA969";
                  };
                })
                .attr("d", function(d){
                  if(d.Strand == 1){
                    return lineFunction(myPathForward(d.Start, d.Stop));
                  }else{
                    return lineFunction(myPathBackward(d.Start, d.Stop));
                  };
                });

  }
