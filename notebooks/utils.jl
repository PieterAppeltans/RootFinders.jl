using PyPlot

function plotf(x,f)
	y = f(x);
	plot(x,y,linewidth=2.);
	xlabel("x",x=1);
	ylabel("f(x)",rotation=0,y=1);
	ax = gca();
	ax.spines["top"].set_visible(false) # Hide the top edge of the axis
	ax.spines["right"].set_visible(false) # Hide the right edge of the axis
	ax.spines["left"].set_position("center") # Move the right axis to the center
	ax.spines["bottom"].set_position("center") # Most the bottom axis to the center
	ax.xaxis.set_ticks_position("bottom") # Set the x-ticks to only the bottom
	ax.yaxis.set_ticks_position("left") # Set the y-ticks to only the left

	
end
