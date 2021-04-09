function stackedBar(aID,dataSet1){
	var width       = $('#STACBARwidth').val(),
	    height      = $('#STACBARheight').val(),
        svgHeight   = Number(height)+100;
	    padding     = {left: 50, right: 200, top:20, bottom: 20},
	    xRangeWidth = width - padding.left - padding.right,
	    yRangeHeight = height - padding.top - padding.bottom;
	
	var vis = d3v3.select("#"+aID).append("div").attr({
	        margin: "auto",
	        id: "STACBARplot"
	    }),
	    svg = vis
	    .append("svg")
	    .attr("id", "STACBARsvg")
        .attr("xmlns", "http://www.w3.org/2000/svg")
        .attr("width", width)
        .attr("height", svgHeight)
	    .append("g")
	    .attr("transform", "translate(" + [padding.left, padding.top] + ")");
	
	var offsetSelect = d3v3.ui.select({
	        base: vis,
	        before: "svg",
	        style: {position: "absolute", left: 0 + "px", top: -40 + "px"},
	        onchange: function() {
	            update(dataSet1)
	        },
	        data: ["number", "proportion", "silhouette"]
	    }),
	    orderSelect  = d3v3.ui.select({
	        base: vis,
	        before: "svg",
	        style: {position: "absolute", left: 100 + "px", top: -40 + "px"},
	        onchange: function() {
	            update(dataSet1)
	        },
	        data: ["inside-out", "default", "reverse"]
	    }),
	    stack        = d3v3.layout.stack()
	        .values(function(d) { return d.sales; })
	        .x(function(d) { return d.year; })
	        .y(function(d) { return d.profit; })
	        .out(function out(d, y0, y) {
	            d.p0 = y0;
	            d.y = y;
	        }
	    );
	
	// x Axis
	var xPadding = {inner: 0.1, outer: 0.3},
	    xScale   = d3v3.scale.ordinal()
	        .rangeBands([0, xRangeWidth], xPadding.inner, xPadding.outer),
	    xAxis    = d3v3.cbPlot.d3Axis()
	        .scale(xScale)
	        .orient("bottom"),
	    gX       = svg.append("g")
	        .attr("class", "x axis")
	        .attr("transform", "translate(0," + yRangeHeight + ")");
	
	// y Axis
	var yAxisScale = d3v3.scale.linear()
	        .range([yRangeHeight, 0]),
	    yAxis      = d3v3.cbPlot.d3Axis()
	        .scale(yAxisScale)
	        .orient("left")
	        .tickSubdivide(2),
	    gY         = svg.append("g")
	        .attr("class", "y axis")
	        .style({"pointer-events": "none", "font-size": "12px"}),
	    yAxisTransition = 1000;
	
	var yPlotScale = d3v3.scale.linear()
	    .range([0, yRangeHeight]);
	
	var color = d3v3.scale.category10();
	
	function update(dataSet) {
	    // create an array of normalised layers and
	    // add the normalised values onto the data
	    var normData     = stack.offset("proportion")(dataSet)
	            .map(stack.values())
	            .map(function(s) {
	                return s.map(function(p) {return p.yNorm = p.y})
	            }),
	        stackedData  = stack.offset(offsetSelect.value())
	            .order(orderSelect.value())(dataSet),
	        maxY         = d3v3.max(stackedData, function(d) {
	            return d3v3.max(d.sales, function(s) {
	                return s.profit + s.p0
	            })
	        }),
	        years        = stackedData[0].sales.map(stack.x()),
	        yearlyTotals = years.reduce(function(t, y) {
	            return (t[y] = d3v3.sum(stackedData, function(o) {
	                return o.sales.filter(function(s) {
	                    return s.year == y
	                })[0].profit
	            }), t)
	        }, {});
	
	    xScale.domain(years);
	    yAxisScale.reset = function(){
	        this.domain([0, offsetSelect.value() == "proportion" ? 1 : maxY])
	            .range([yRangeHeight, 0])
	            .ticks(10)
	    };
	    yAxisScale.reset();
	    yPlotScale.domain(yAxisScale.domain());
	
	    // plotArea
	    // (svg) -> (g.plotArea)[stackedData]
	    // apply a transform to map screen space to cartesian space
	    // this removes all confusion and mess when plotting data!
	    var plotArea = svg.selectAll(".plotArea")
	        .data([stackedData]);
	    plotArea.enter().insert("g", ".axis")
	        .attr(d3v3.cbPlot.transplot(yRangeHeight))
	        .attr("class", "plotArea");
	
	    plotArea.series = plotArea.selectAll(".series")
	        .data(ID);
	    plotArea.series.enter()
	        .append("g")
	        .attr("class", "series");
	    plotArea.series.style("fill", function(d, i) {
	        return color(i);
	    });
	    plotArea.series.exit().remove();
	    Object.defineProperties(plotArea.series, d3v3._CB_selection_destructure);
	
	    plotArea.series.components = plotArea.series.selectAll(".components")
	        .data(function(d) {
	            return d3v3.entries(d);
	        });
	    plotArea.series.components.enter().append("g")
	        .attr("class", function(d){return d.key})
	        .classed("components", true);
	    plotArea.series.components.exit().remove();
	
	    plotArea.series.components.values = plotArea.series.components.filter(function(d){
	        return d.key == "sales"
	    });
	    Object.defineProperties(plotArea.series.components.values, d3v3._CB_selection_destructure);
	    plotArea.series.components.labels = plotArea.series.components.filter(function(d){
	        return d.key == "name"
	    })
	        // reverse the plotArea transform (it is it's own inverse)
	        .attr(d3v3.cbPlot.transplot(yRangeHeight));
	    Object.defineProperties(plotArea.series.components.labels, d3v3._CB_selection_destructure);
	
	    var s         = xScale.rangeBand(),
	        w         = s - xPadding.inner,
	        drag = d3v3.behavior.drag()
	            .on("dragstart", mouseOver),
	        points = plotArea.series.components.values.points = plotArea.series.components.values.selectAll("rect")
	            .data(function(d){
	                return d.value
	            });
	    points.enter()
	        .append("rect")
	        .attr({width: w, class: "point"})
	        .on("mouseover", mouseOver)
	        .on("mouseout", mouseOut)
	        .call(drag);
	    points.transition()
	        .attr("x", function(d) {
	            return xScale(d.year);
	        })
	        .attr("y", function(d) {
	            return yPlotScale(d.p0);
	        })
	        .attr("height", function(d) {
	            return yPlotScale(d.y);
	        })
	        .attr("stroke", "white");
	
	    points.exit().remove();
	    Object.defineProperties(plotArea.series.components.values.points, d3v3._CB_selection_destructure);
	
		svg.selectAll(".x.axis .tick text").style("text-anchor", "end").attr("dx",$('#STACBARxlabelshift').val()+"px").attr("transform", "rotate("+$('#STACBARxlabelrotate').val()+")").style("font-size",$('#STACBARxfontsize').val()+"px");;
	
	    gX.transition().call(xAxis);
	    gY.transition().call(yAxis);
	
	    function mouseOver(pointData, pointIndex, groupIndex) {
// console.log(["in", pointIndex].join("\t"));
	        var selectedYear = pointData.year,
	            // wrap the node in a selection with the proper parent
	            plotData = plotArea.series.components.values.data,
	            seriesData = plotData[groupIndex],
	            currentYear = d3v3.transpose(plotData)[pointIndex],
	            point         = plotArea.series.components.values.points.nodes[groupIndex][pointIndex];
	
	        // if the plot is not normalised, fly-in the axis on the selected year
	        if(offsetSelect.value() != "proportion") {
	            yAxisScale.reset();
	            // get the zero offset for the fly-in axis
	            var pMin        = d3v3.min(currentYear, function(s) {
	                    return s.p0
	                }),
	                refP0 = seriesData[pointIndex].p0,
	                selectedGroupHeight = d3v3.sum(currentYear, function(d) {return d.y}),
	                // set the range and domain height for the selected year
	                localDomain = [0, selectedGroupHeight].map(function(d){return d + pMin - refP0}),
	                localRange  = [0, selectedGroupHeight].map(function(d) {return yAxisScale(d + pMin)});
// console.log(yAxisScale(pMin));
	            yAxisScale
	                .domain(localDomain)
	                .range(localRange);
	            // apply the changes to the y axis and manage the ticks
	            gY.transition("axis")
	                .duration(yAxisTransition)
	                .call(yAxis.ticks(+(Math.abs(localRange[0] - localRange[1]) / 15).toFixed()))
	                .attr("transform", "translate(" + point.attr("x") + ",0)")
	                .style({"font-size": "8px"})
	                .call(function(t) {d3v3.select(t.node()).classed("fly-in", true)});
	            // align the selected series across all years
	            points.transition("points")
	                .attr("y", alignY(seriesData[pointIndex].p0, groupIndex))
	                .call(endAll, toolTip)
	        } else window.setTimeout(toolTip, 0);  // if not proportion
	
	        // manage the highlighting
	        //  points highlighting
	        plotArea.series.transition("fade")
	            .attr("opacity", function(d, i) {
	                return i == groupIndex ? 1 : 0.5;
	            });
	        //  x axis highlighting
	        d3v3.selectAll(".x.axis .tick")
	            .filter(function(d) {
	                return d == selectedYear
	            })
	            .classed("highlight", true);
	
	        // move the selected element to the front
	        d3v3.select(this.parentNode)
	            .moveToFront();
	        gX.moveToFront();
	
	
	        legendText(groupIndex);
	
	        // Tooltip
	        function toolTip() {
	            plotArea.series
	                .append("g")
	                .attr("class", "tooltip")
	                .attr("transform", "translate(" + [point.attr("x"), point.attr("y")] + ")")
	                .append("text")
	                .attr(d3v3.cbPlot.transflip())
	                .text(d3v3.format(">8.0%")(pointData.yNorm))
	                .attr({x: "1em", y: -point.attr("height") / 2, dy: ".35em", opacity: 0})
	                .transition("tooltip").attr("opacity", 1)
	                .style({fill: "black", "pointer-events": "none"})
	        }
	    }
	    function mouseOut(d, nodeIndex, groupIndex) {
// console.log(["out", nodeIndex].join("\t"));
	        var year = d.year;
	        d3v3.selectAll(".x.axis .tick")
	            .filter(function(d) {
	                return d == year
	            })
	            .classed("highlight", false);
	        plotArea.series.transition("fade")
	            .attr({opacity: 1});
	        var g = plotArea.series.components.labels.nodes[groupIndex][0].select("text");
	        g.classed("highlight", false);
	        g.text(g.text().split(":")[0])
	        yAxisScale.reset();
	        gY.selectAll(".minor").remove();
	        gY.transition("axis").call(yAxis)
	            .attr("transform", "translate(0,0)")
	            .style({"font-size": "12px"})
	            .call(function(t) {d3v3.select(t.node()).classed("fly-in", false)});
	        plotArea.series.selectAll(".tooltip")
	            .transition("tooltip")
	            .attr({opacity: 0})
	            .remove();
	        points.transition("points").attr("y", function(d) {
	            return yPlotScale(d.p0);
	        })
	    };
	
	
	    // Add the legend inside the series containers
	    // The series legend is wrapped in another g so that the
	    // plot transform can be reversed. Otherwise the text would be mirrored
	    var labHeight = 20,
	        labRadius = 6;
	
	    // add the marker and the legend text to the normalised container
	    // push the stackedData (name) down to them
	    var labelCircle = plotArea.series.components.labels.selectAll("circle")
	            .data(function(d){return [d.value]}),
	        // take a moment to get the series order delivered by stack
	        orders      = stackedData.map(function(d) { // simplify the form
	            return {name: d.name, base: d.sales[0].p0}
	        }).sort(function(a, b) {        // get a copy, sorted by p0
	            return a.base - b.base
	        }).map(function(d) {            // convert to index permutations
	            return stackedData.map(function(p) {
	                return p.name
	            }).indexOf(d.name)
	        }).reverse();                   // convert to screen y ordinate
	    labelCircle.enter().append("circle")
	        .on("mouseover", function(pointData, pointIndex, groupIndex) {
	            var node = this,
	                typicalP0 = d3v3.median(plotArea.series.components.values.data[groupIndex],
	                function(d){return d.p0});
	            plotArea.series.components.values.points.transition("points")
	                .attr("y", alignY(typicalP0, groupIndex));
	            plotArea.series.transition("fade")
	                .attr("opacity", function(d) {
	                    return d === d3v3.select(node.parentNode.parentNode).datum() ? 1 : 0.5;
	                });
	            legendText(groupIndex);
	        })
	        .on("mouseout", function(pointData, pointIndex, groupIndex) {
	            plotArea.series.transition("fade")
	                .attr({opacity: 1});
	            plotArea.series.components.values.points.transition("points").attr("y", function(d) {
	                return yPlotScale(d.p0);
	            })
	        });
	    labelCircle.attr("cx", xRangeWidth + 10)
	        .attr("cy", function(d, i, j) {
	            return labHeight * orders[j];
	        })
	        .attr("r", labRadius);
	
	    var labelText = plotArea.series.components.labels.selectAll("text")
	        .data(function(d){return [d.value]});
	    labelText.enter().append("text");
	    labelText.attr("x", xRangeWidth + 20)
	        .attr("y", function(d, i, j) {
	            return labHeight * orders[j];
	        })
	        .attr("dy", labRadius / 2)
	        .text(function(d) {
	            return d;
	        });
	
	    function legendText(groupIndex){
	        // Legend text
	        // add the value for the moused over item to the legend text and
	        // highlight it
	        var labelText = plotArea.series.components.labels.nodes[groupIndex][0].select("text"),
	            seriesData = plotArea.series.components.values.data[groupIndex],
	            fmt           = [">8,.0f", ">8.0%"][(offsetSelect.value() == "proportion") * 1];
	        labelText.classed("highlight", true);
	        labelText.text(labelText.datum().value + ": " + d3v3.format(fmt)(
	                offsetSelect.value() != "proportion" ?
	                d3v3.sum(seriesData, stack.y()) :
	                d3v3.sum(seriesData, function(s) {
	                    var totalSales = d3v3.sum(d3v3.values(yearlyTotals));
	                    return s.y * yearlyTotals[s.year] / totalSales
	                })
	            ));
	    }
	    function alignY(p0, series) {
	        var offsets = plotArea.series.components.values.data[series].map(function(d) {
	            return p0 - d.p0;
	        });
	        return function(d, i) {
	            return yPlotScale(d.p0 + offsets[i]);
	        }
	    }
	
	    function aID(d) {
	        return [d];
	    }
	    function ID(d) {
	        return d;
	    }
	}
	
	d3v3.selection.prototype.moveToFront = function() {
	    return this.each(function() {
	        this.parentNode.appendChild(this);
	    });
	};
	d3v3._CB_selection_destructure = {
	    "nodes": {
	        get: function() {
	            return this.map(function(g) {
	                return g.map(function(n) {
	                    return d3v3.select(n)
	                })
	            })
	        }
	    },
	    data: {
	        get: function() {
	            return this.map(function(g) {
	                return d3v3.select(g[0]).datum().value
	            })
	        }
	    }
	};
	
//    svg.append("style").text("path { fill: none; stroke: #000; shape-rendering: crispEdges; }")
	update(dataSet1);
	svg.selectAll(".axis path").style("fill","none").style("stroke","#000").style("shape-rendering","crispEdges");

	window.setTimeout(function(){
	    update(dataSet1.map(function(d) {
	        return {
	            name: d.name, sales: d.sales.map(function(y) {
	                return {year: y.year, profit: y.profit / 2}
	            })
	        }
	    })
	    )
	},1000)
}
