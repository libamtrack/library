package amtrackservice.client.gui.elements;

import java.util.HashMap;

import amtrackservice.client.MapList;
import ca.nanometrics.gflot.client.DataPoint;
import ca.nanometrics.gflot.client.PlotItem;
import ca.nanometrics.gflot.client.PlotModel;
import ca.nanometrics.gflot.client.PlotPosition;
import ca.nanometrics.gflot.client.SeriesHandler;
import ca.nanometrics.gflot.client.SimplePlot;
import ca.nanometrics.gflot.client.event.PlotHoverListener;
import ca.nanometrics.gflot.client.jsni.Plot;
import ca.nanometrics.gflot.client.options.AxisOptions;
import ca.nanometrics.gflot.client.options.GridOptions;
import ca.nanometrics.gflot.client.options.LineSeriesOptions;
import ca.nanometrics.gflot.client.options.PlotOptions;
import ca.nanometrics.gflot.client.options.PointsSeriesOptions;
import com.google.gwt.user.client.ui.FlowPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.Label;
import com.google.gwt.user.client.ui.Widget;

public class AmPlot extends AmWidget {

	private static final String INSTRUCTIONS = "Point your mouse to a data point on the chart";

	private final Label selectedPointLabel = new Label(INSTRUCTIONS);

	private MapList<String, String> preset;
	private String dataX;
	private String dataY;
	private HTML description;
	private FlowPanel widget = new FlowPanel();

	public AmPlot(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX, String dataY,
			String dataZ) {
		super(label, datatype, description);
		this.preset = preset;
		this.dataX = dataX;
		this.dataY = dataY;
		widget.add(new HTML("No data"));

	}

	@Override
	public String getValue() {
		// unimplemented
		return "";
	}

	@Override
	public Widget getWidget() {
		return widget;
	}

	@Override
	public void setValue(HashMap<String, String> valueMap) {
		HashMap<Double, Double> values = new HashMap<Double, Double>();
		String[] xValues = valueMap.get(dataX).split(" ");
		String[] yValues = valueMap.get(dataY).split(" ");
		for (int i = 0; (i < xValues.length) && (i < yValues.length); i++) {
			double x = Double.parseDouble(xValues[i]);
			double y = Double.parseDouble(yValues[i]);
			values.put(x, y);
		}
		this.widget.clear();
		this.widget.add(createPlot(values));
	}

	@Override
	public String getDataLink() {
		return null; // unimplemented
	}

	private Widget createPlot(final HashMap<Double, Double> values) {
		PlotModel model = new PlotModel();
        PlotOptions plotOptions = new PlotOptions();
        plotOptions.setDefaultLineSeriesOptions(new LineSeriesOptions().setLineWidth(1).setShow(true));
        plotOptions.setDefaultPointsOptions(new PointsSeriesOptions().setRadius(2).setShow(true));
        plotOptions.setDefaultShadowSize(0);

        plotOptions.setXAxisOptions(new AxisOptions() {
		
        });
        
        // >>>>>>> You need make the grid hoverable <<<<<<<<<
        plotOptions.setGridOptions(new GridOptions().setHoverable(true));
        // >>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        // create a series
        SeriesHandler handler = model.addSeries("Data Line", "#007f00");

        for(Double d: values.keySet()){
        	handler.add(new DataPoint(d,values.get(d)));
        }

        // create the plot
        SimplePlot plot = new SimplePlot(model, plotOptions);

        // add hover listener
        plot.addHoverListener(new PlotHoverListener() {
                public void onPlotHover(Plot plot, PlotPosition position, PlotItem item) {
                        if (item != null) {
                                selectedPointLabel
                                                .setText("x: " + item.getDataPoint().getX() + ", y: " + item.getDataPoint().getY());
                        } else {
                                selectedPointLabel.setText(INSTRUCTIONS);
                        }
                }
        }, false);

        // put it on a panel
        FlowPanel panel = new FlowPanel();
        panel.add(plot);
        panel.add(selectedPointLabel);

        return panel;

	}

	@Override
	public void setDefault() {
		// TODO Auto-generated method stub
		
	}

}
