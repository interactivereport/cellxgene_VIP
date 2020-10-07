function volcanoPlot() {
  var width = 960,
    height = 500,
    margin = {top: 20, right: 20, bottom: 50, left: 50},
    xColumn, // name of the variable to be plotted on the axis
    yColumn,
    sizeColumn="",
    colorColumn="",
    labelColumn="",
    xAxisLabel, // label for the axis
    yAxisLabel,
    sampleID = "Gene",
    hLines={},
    vLines={},
    xScale = d3.scaleLinear(), // the values for the axes will be continuous
    yScale = d3.scaleLinear(),
    maxY = 20,
    maxX = 10,
    tooltips=[],
    resetBtnID = ""; //the reset button element id
    colGrp={};
    
    function chart(selection){
      var innerWidth = width - margin.left - margin.right, // set the size of the chart within its container
          innerHeight = height - margin.top - margin.bottom;

      selection.each(function(data) {
        //check the color
        if(colorColumn.length>1){
          var colGrpName = d3.map(data,function(d){return d[colorColumn];}).keys();
          colGrp = {};
          if(colGrpName.includes('unselected')){
            colGrp['unselected']='dot Unselect';
            colGrpName=colGrpName.filter(function(e) { return e !== 'unselected' });
          }
          var n=colGrpName.length;
          while(n--) colGrp[colGrpName[n]] = 'dot col'+(n+1);
        }

            // set up the scaling for the axes based on the inner width/height of the chart and also the range
            // of value for the x and y axis variables. This range is defined by their min and max values as
            // calculated by d3.extent()
            var xlim = d3.extent(data, function(d) { return d[xColumn]; });
            xlim[0] -= 0.1*(xlim[1]-xlim[0]);
            xlim[1] += 0.1*(xlim[1]-xlim[0]);
            xScale.range([0, innerWidth])
                .domain(xlim)
                .nice();

            // normally would set the y-range to [height, 0] but by swapping it I can flip the axis and thus
            // have -log10 scale without having to do extra parsing
            var ylim = d3.extent(data, function(d) { return yNegLog(d[yColumn]); });
            ylim[0] -= 0.1*(ylim[1]-ylim[0]);
            ylim[1] += 0.1*(ylim[1]-ylim[0]);
            yScale.range([innerHeight,0])
                .domain(ylim)
                .nice(); // adds "padding" so the domain extent is exactly the min and max values

            var zoom = d3.zoom()
                .scaleExtent([1, 20])
                .translateExtent([[0, 0], [width, height]])
                .on('zoom', zoomFunction);

            // append the svg object to the selection
            var svg = d3.select(this).append('svg')
                .attr('height', height)
                .attr('width', width)
                .append('g')
                .attr('transform', 'translate(' + margin.left + ',' + margin.top + ')')
                .call(zoom);

            // position the reset button and attach reset function
            if(resetBtnID.length>2){
              d3.select('#'+resetBtnID)
                .attr('x', margin.top * 1.5 + 'px')
                .attr('y', margin.left * 1.25 + 'px')
                .on('click', reset);
            }

            svg.append('defs').append('clipPath')
                .attr('id', 'clip')
                .append('rect')
                .attr('height', innerHeight)
                .attr('width', innerWidth);

            // add the axes 
            var xAxis = d3.axisBottom(xScale);
            var yAxis = d3.axisLeft(yScale)
                          .ticks(5);
                          //.tickFormat(yTickFormat);
            var gX = svg.append('g')
                .attr('class', 'x axis')
                .attr('transform', 'translate(0,' + innerHeight + ')')
                .call(xAxis);
            gX.append('text')
                .attr('class', 'label')
                .attr('transform', 'translate(' + innerWidth / 2 + ',' + (margin.bottom - 12) + ')')
                .style('text-anchor', 'middle')
                .html(xAxisLabel || xColumn);
            var gY = svg.append('g')
                .attr('class', 'y axis')
                .call(yAxis);
            gY.append('text')
                .attr('class', 'label')
                .attr('transform', 'translate(' + (0 - margin.left / 1.25) + ',' + (height / 2) + ') rotate(-90)')
                .style('text-anchor', 'middle')
                .html(yAxisLabel || yColumn);

            // this rect acts as a layer so that zooming works anywhere in the svg. otherwise, if zoom is called on
            // just svg, zoom functionality will only work when the pointer is over a circle.
            var zoomBox = svg.append('rect')
                .attr('class', 'zoom')
                .attr('height', innerHeight)
                .attr('width', innerWidth);

            var circles = svg.append('g')
                .attr('class', 'circlesContainer');

            circles.selectAll(".dot")
              .data(data)
              .enter().append('circle')
              .attr("r", circleSize)
              .attr('cx', function(d) { return xScale(d[xColumn]); })
              .attr('cy', function(d) { return yScale(yNegLog(d[yColumn])); })
              .attr('class', circleClass)
              .on('mouseenter', tipEnter)
              .on("mousemove", tipMove)
              .on('mouseleave', function(d) {
                 return tooltip.style('visibility', 'hidden');
              });
              
            
            //add labels
            var dragLabel = d3.drag()
              .on("drag", function(d) {
                circles.select("#point"+d[sampleID])
                  .attr("x2",transform.invertX(d3.event.x))//
                  .attr("y2",transform.invertY(d3.event.y));
                d3.select(this)
                  .attr("x", transform.invertX(d3.event.x))
                  .attr("y", transform.invertY(d3.event.y))
            });
            circles.selectAll("text")
              .data(data.filter(function(d){return d[labelColumn];}))
              .enter()
              .append('text')
              .attr('x', function(d) { return xScale(d[xColumn]); })
              .attr('y', function(d) { return yScale(yNegLog(d[yColumn])); })
              .attr("class", "dotName")
              .text(function(d){return d[sampleID];})
              .call(dragLabel)
            circles.selectAll("line")
              .data(data.filter(function(d){return d[labelColumn];}))
              .enter()
              .append("line")
              .attr("x1", function(d) { return xScale(d[xColumn]); })
              .attr("x2", function(d) { return xScale(d[xColumn]); })
              .attr("y1", function(d) { return yScale(yNegLog(d[yColumn])); })
              .attr("y2", function(d) { return yScale(yNegLog(d[yColumn])); })
              .attr("id",function(d){return 'point'+d[sampleID]})
              .attr('class','dotLine');

            /* adding lines */
            var thresholdLines = svg.append('g')
                .attr('class', 'thresholdLines');
            // add horizontal lines
            for (const [key, value] of Object.entries(hLines)) {
              thresholdLines.append("svg:line")
                  .attr('class', 'threshold')
                  .attr("x1", 0)
                  .attr("x2", innerWidth)
                  .attr("y1", yScale(yNegLog(value)))
                  .attr("y2", yScale(yNegLog(value)));
              if(key.length>1){
                thresholdLines.append("text")
                    .attr('class','zoomLabel')
                    .attr("x", xScale(0))
                    .attr("y", yScale(yNegLog(value)))//*1.15
                    .attr("fill", "red")
                    .style('dominant-baseline', 'hanging')
                    .style('text-anchor', 'middle')
                    .text(key);
              }                
            }

            // add vertical lines
            for (const [key, value] of Object.entries(vLines)) {
              //console.log(key, value);
              thresholdLines.append("svg:line")
                  .attr('class', 'threshold')
                  .attr("x1", xScale(value))
                  .attr("x2", xScale(value))
                  .attr("y1", 1)
                  .attr("y2", innerHeight);
              if(key.length>1){
                thresholdLines.append("text")
                    .attr('class','zoomLabel')
                    .attr("x", xScale(value))
                    .attr("y", 0)
                    .attr("fill", "red")
                    .style('text-anchor', 'middle')
                    .text(key);
              }                
            }
            
            //add tooltip
            var tooltip = d3.select("body")
                .append("div")
                .attr('class', 'volcanoTooltip');

            function tipEnter(d) {
              strHtml = "";
              tooltips.forEach(function(one){
                strHtml += '<strong>'+one+'</strong>: '+d[one]+'<br/>'
              })
              //console.log(strHtml);
              tooltip.style('visibility', 'visible')
                      .style('font-size', '11px')
                      .html(strHtml);
            }

            function tipMove() {
                tooltip.style("top", (event.pageY - 5) + "px")
                    .style("left", (event.pageX + 20) + "px");
            }

            function yTickFormat(n) {
                return d3.format(".2r")(getBaseLogNeg(10, n));
                function getBaseLog(x, y) {
                    return Math.log(y) / Math.log(x);
                }
                function getBaseLogNeg(x, y) {
                    return -Math.log(y) / Math.log(x);
                }
            }
            function yNegLog(y){
              var tmp = Math.round(100*(Number.EPSILON-Math.log(y)/Math.log(10)))/100;
              //if(!isFinite(tmp))console.log(isFinite(tmp) ? Math.min(tmp,maxY):maxY);
              return isFinite(tmp) ? Math.min(tmp,maxY):maxY;
            }

            var transform = d3.zoomIdentity;
            function zoomFunction() {
                transform = d3.zoomTransform(this);
                d3.selectAll('.dot')
                    .attr('transform', transform)
                    .attr('r', circleSize);//Math.sqrt(2 / transform.k
                gX.call(xAxis.scale(d3.event.transform.rescaleX(xScale)));
                gY.call(yAxis.scale(d3.event.transform.rescaleY(yScale)));
                svg.selectAll('.threshold')
                    .attr('transform', transform)
                    .attr('stroke-width', 1 / transform.k);
                    
                svg.selectAll('.zoomLabel')
                    .attr('transform', transform);
                svg.selectAll('.dotName')
                    .attr('transform', transform)
                    .attr('font-size',16/transform.k);
                svg.selectAll('.dotLine')
                    .attr('transform', transform);
            }

            function circleClass(d) {
              if(colorColumn.length<2) return 'dotUnselect';
              return colGrp[d[colorColumn]];
            }
            function circleSize(d) {
              var scale = 1;
              if(typeof transform != "undefined"){
                scale = transform.k
              }
              if(sizeColumn.length<2) return 3.5/scale;
              if (d[sizeColumn] < 3.5) {
                return 3.5/scale;
              } else {
                return d[sizeColumn]/scale;
              }
            }

            function reset() {
              var ease = d3.easePolyIn.exponent(4.0);
              svg.transition().duration(750)
                .ease(ease)
                .call(zoom.transform, d3.zoomIdentity);
            }
        });
    }

    chart.width = function(value) {
        if (!arguments.length) return width;
        width = value;
        return chart;
    };

    chart.height = function(value) {
        if (!arguments.length) return height;
        height = value;
        return chart;
    };

    chart.margin = function(value) {
        if (!arguments.length) return margin;
        margin = value;
        return chart;
    };

    chart.xColumn = function(value) {
        if (!arguments.length) return xColumn;
        xColumn = value;
        return chart;
    };

    chart.yColumn = function(value) {
        if (!arguments.length) return yColumn;
        yColumn = value;
        return chart;
    };
    chart.sizeColumn = function(value) {
        if (!arguments.length) return sizeColumn;
        sizeColumn = value;
        return chart;
    };
    chart.colorColumn = function(value) {
        if (!arguments.length) return colorColumn;
        colorColumn = value;
        return chart;
    };
    chart.labelColumn = function(value) {
        if (!arguments.length) return labelColumn;
        labelColumn = value;
        return chart;
    };
    
    chart.xAxisLabel = function(value) {
        if (!arguments.length) return xAxisLabel;
        xAxisLabel = value;
        return chart;
    };

    chart.yAxisLabel = function(value) {
        if (!arguments.length) return yAxisLabel;
        yAxisLabel = value;
        return chart;
    };

    chart.sampleID = function(value) {
        if (!arguments.length) return sampleID;
        sampleID = value;
        return chart;
    };
    
    chart.resetBtnID = function(value) {
        if (!arguments.length) return resetBtnID;
        resetBtnID = value;
        return chart;
    };
    chart.maxY = function(value) {
        if (!arguments.length) return maxY;
        maxY = value;
        return chart;
    };
    chart.vLines = function(value) {
        if (!arguments.length) return vLines;
        vLines = value;
        return chart;
    };
    chart.hLines = function(value) {
        if (!arguments.length) return hLines;
        hLines = value;
        return chart;
    };
    chart.tooltips = function(value) {
        if (!arguments.length) return tooltips;
        tooltips = value;
        return chart;
    };
    return chart;
}