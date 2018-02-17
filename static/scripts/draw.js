function draw(parsedData){
  var margin = {top: 20, right: 100, bottom: 20, left: 100};

  width = document.documentElement.clientWidth;

  var width = width - margin.left - margin.right;
  var height = 200 - margin.top - margin.bottom;

  var svgContainer = d3.select(".svg").append("svg")
                        .attr("width", width + margin.left + margin.right)
                        .attr("height", height + margin.top + margin.bottom)
                        .append("g")
                        .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
                        .attr("style", "outline: thin solid red")
                        .attr("viewbox", "0 0 1500 200")
                        .attr("preserveAspectRatio","xMidYMid meet");

  parsedData.sort(function(a, b){
    return a.Start - b.Start
  });

  // Define the x-axis horizontal
  var minValue = parsedData[0].Start;
  var maxValue = 5000 + parsedData[parsedData.length - 1].Stop;

  xScale = d3.scaleLinear()
                       .domain([minValue, maxValue])
                       .range([100, width]);

  var line = svgContainer.selectAll("line")
                .data(parsedData)
                .enter().append("svg:line")
                .attr("x1", function(d, i){
                  if(i== 0){
                    return xScale(d.Start)
                  }
                  return xScale(d.Start)})
                .attr("y1", "50")
                .attr("x2", function(d){ return xScale(d.Stop)})
                .attr("y2", "50")
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

  }

  d3.select(window).on('resize', resize);

  function resize(){
    // update width
    width = document.documentElement.clientWidth;

    // reset x range

  }
