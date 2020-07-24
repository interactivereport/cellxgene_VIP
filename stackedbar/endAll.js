function endAll(transition, callback) {
	var n = 0;
	transition
		.each(function() { ++n; })
		.each("end.endAll", function() { if (!--n) callback.apply(this, arguments); });
	if(transition.empty()) callback();
}
