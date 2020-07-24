(function(d3) {

    var exports = d3.cbPlot = d3.coolPlot || {};

    /**
     * Created by cool.blue@y7mail.com on 22/08/2015.
     * Returns an afine transform that maps from top left to bottom left origin
     *  the transform is wrapped in an object that will be accepted by .attr() in d3
     * @param height
     *    - height of the plot area in pixels
     * @returns
     *  - attr object
     */

    function transplot(yRange) {
        return {"transform": "matrix(" + [1, 0, 0, -1, 0, yRange] + ")"};
    }

    /**
     * Reverses the local mirroring of transplot
     * @returns {{transform: string}}
     */
    function transflip() {
        return {"transform": "matrix(" + [1, 0, 0, -1, 0, 0] + ")"};
    }

    /**
     *
     * @param tickSize
     * @returns {tickSize}
     */
    function tickSize(tickSize) {
        var axis      = this,
            tickSize0 = Math.max(axis.innerTickSize(), 0),
            tickSize1 = Math.max(tickSize, 0),
            padding   = Math.max(axis.tickPadding(), 0) + tickSize0 - tickSize1;
        axis.innerTickSize(tickSize).tickPadding(padding);
        return this;
    }

    function tickSubdivide(n) {
        var axis = this, scale = axis.scale(),
            tickLabels = scale.ticks,
            ε = 1e-6;

        function coolAxis(selection) {
            var minorTickValues = [], scale0;
            selection.each(function(){
                scale0 = this.__chart__ || axis.scale();
                d3.select(this).selectAll(".tick.minor").each(function(d, i){
                    minorTickValues.push(d);
                }).remove()
            });
            console.log(minorTickValues.map(f).join(""))
            selection.call(axis);
            // use each and reselect to allow for transition objects as well as selections
            selection.each(function() {
                var g = d3.select(this);
                g.selectAll(".tick.minor").data(minorTickValues)
                    .enter().insert("g", ".domain")
                    .attr({
                        class: "tick minor",
                        transform: function(d) {return "translate(0," + scale0(d) + ")"}
                    }).append("line")
                    .attr({
                        "y2": 0,
                        "x2": -3
                    });

                var ticks = g.selectAll(".tick")
                        .data(scale.copy().ticks(axis.ticks()[0] * n), scale),
                    enter = ticks.enter().insert("g", ".domain")
                        .attr({
                            class: "tick minor"
                        })
                        .style("opacity", ε),
                    exit = d3.transition(ticks.exit())
                        .style("opacity", ε).remove(),
                    update = d3.transition(ticks.order())
                        .style("opacity", 1);
                d3.transition(enter)
                    .attr({
                        transform: function(d) {return "translate(0," + scale(d) + ")"}
                    })
                enter.append("line")
                    .transition("minor")
                    .attr({
                        "y2": 0,
                        "x2": -3
                    })
            });

            function tlog(trans, type, message){
                var data = Array(trans[0].length).join("_").split("");
                console.log(padRight(message, 12) + (trans.each(function(d, i){
                    data[i] = d;
                }), data.map(f).join(" ")));
                if(trans.duration && type) {
                    type = Array.isArray(type) ? type : [type];
                    type.forEach(function(t) {
                        trans.each(t, function(d, i, j) {
                            var tag = d3.select(this);
                            console.log([t, "on", tag.attr("class") || tag.property("tagName"), "x", trans.size()].join(" "))
                        })
                    })
                }
            }
            function f(x){
                var l = 8;
                return x == null || x == undefined || x == "_" ? Array(l+1).join(".") : d3.format("_>"+l+"d")(x);
            }
            function padRight(s, l){
                l = l || 8;
                var len = s.length;
                return s + Array(l - len).join("_");
            }
        }

        return d3.rebind.bind(null, coolAxis, axis).apply(null, Object.keys(axis));
    }

    /**
     * Axis constructor that returns custom behaviour on d3.svg.axis
     *
     */
    function d3TransfAxis() {
        var axis = d3.svg.axis();

        function transAxis(g) {
            g.call(axis).selectAll(".tick text, .tick line").attr(transflip())
            if(d3.select(g.node()).classed("x")) g.selectAll(".domain").attr(transflip());
        }

        axis.tickSize = tickSize.bind(axis);
        d3.rebind.bind(null, transAxis, axis).apply(null, Object.keys(axis))
        return transAxis;
    }
    exports.transplot = transplot;
    exports.transflip = transflip;

    exports.d3Axis = function d3Axis() {
        var axis = d3.svg.axis();
        axis.tickSize = tickSize.bind(axis);
        axis.tickSubdivide = tickSubdivide.bind(axis);
        return axis;
    }
})(d3)