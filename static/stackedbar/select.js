(function(d3) {
	// Example
	// <script src="https://gitcdn.xyz/repo/cool-Blue/d3-lib/master/inputs/select/select.js"></script>
	//
	//		isoLines = d3.ui.select({
	//			base: inputs,
	// 			onUpdate: update,
	//			data: [{text: "show lines", value: "#ccc"}, {text: "hide lines", value: "none"}]
	//		}),
	//
	//		.style("stroke", isoLines.value());
	//

	d3.ui = d3.ui || {};
	d3.ui.select = function (config) {
        // add a select element on base with options matching data
        // if the text and value is the same then data is scalar array
        // 	otherwise the data elements must have text and value fields
        // config
        //  base
        //  before
        //  style
        //  initial
        //  hook
        //  onX
        var select  = (config.base ?
                       (config.base.append ? config.base : d3.select(config.base)) :
                       d3.select("body"))
                [config.before ? "insert" : "append"]("select", config.before ? config.before : null)
                .each(hookEvents)
                .data([config.data]),
            options = select.selectAll("option").data(function(d) {return d});
        options.enter().append("option");
        options.exit().remove();
        if(config.style) select.style(config.style);
        merge(config, options, ["base", "before", "style", "initial", "hook", /on.+/]);
        return options
            .attr({
                value: function(d) {
                    return d.value || d;
                },
                selected: function(d){return d == config.initial ? "selected" : null}
            })
            .text(function(d) {
                return d.text || d
            })
            .call(function() { //add a custom property to the final selection
                if(config.hook) config.hook();
                this.value = function() {
                    return this[0].parentNode.value
                }
            });

        function hookEvents() {
            // store the DOM element
            var _control = this;
            // config object for the control
            // parse the keys and bind a listener to any key beginning with "on"
            Object.keys(config).filter(function(k) {
                return k.slice(0, 2) == "on";
            }).map(function(p, i, listeners) {
                // strip the event name off the listener
                var e = p.slice(2);

                if(!config.on)
                // lazily create a dispatch with the events found in config
                    config.on = d3.dispatch
                        .apply(null, listeners.map(function(l) {
                            return l.slice(2);
                        }));

                config.on.on(e, config[p]);
                // then hook the dispatch to the element event
                d3.select(_control).on(e, function() {
                    config.on[e].apply(this, arguments)
                })

            });
        }
        function merge(source, target, exclude) {
            function included(test){
                return !exclude.some(function(p,i){return test.match(p)})
            }
            for(var p in source)
                if(included(p) && target && !target.hasOwnProperty(p))
                    Object.defineProperty(target, p, Object.getOwnPropertyDescriptor(source, p));
            return target;
        }
    }

    //experiment with chaining...
	d3.ui.select2 = function (base){
		var _base, _on, _onUpdate, _data;
		function select(){}
		d3_ui_select.base = function(base){
			_base = (base.append ? base : d3.select(base)).append("select");
			return this;
		}
		function on(on, listener) {
			if (arguments.length == 1) return _on;
			return this.on(_on || "input", _onUpdate)
		}
	}

})(d3)
