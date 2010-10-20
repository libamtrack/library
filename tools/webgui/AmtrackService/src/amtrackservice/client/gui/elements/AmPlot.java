package amtrackservice.client.gui.elements;

import java.util.HashMap;

import amtrackservice.client.MapList;
import amtrackservice.client.gui.elements.ScatterChartLogAxis.AxisOptions;
import amtrackservice.client.gui.elements.ScatterChartLogAxis.Options;

import com.google.gwt.user.client.ui.FlowPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.AbstractDataTable.ColumnType;
import com.google.gwt.visualization.client.DataTable;

public class AmPlot extends AmWidget {

	private String dataX;
	private String dataY;
	private FlowPanel widget = new FlowPanel();
	
	private HashMap<Double, Double> values = new HashMap<Double, Double>();

	public AmPlot(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX, String dataY,
			String dataZ) {
		super(label, datatype, description);
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
		this.values.clear();
		String[] xValues = valueMap.get(dataX).split(" ");
		String[] yValues = valueMap.get(dataY).split(" ");
		for (int i = 0; (i < xValues.length) && (i < yValues.length); i++) {
			double x = Double.parseDouble(xValues[i]);
			double y = Double.parseDouble(yValues[i]);
			this.values.put(x, y);
		}
		this.widget.clear();
		this.widget.add(createPlot());
	}

	@Override
	public String getDataLink() {
		return null; // unimplemented
	}
		
	private Widget createPlot() {
        FlowPanel panel = new FlowPanel();        
        ScatterChartLogAxis line = new ScatterChartLogAxis(createTable(), createOptions() );
        panel.clear();
        panel.add(line);
        return panel;

	}		

	private Options createOptions() {
		Options options = Options.create();
		options.setWidth(640);
		options.setHeight(480);
		options.setTitle(this.getLabel().getText());
		options.setLineSize(1);
		
		AxisOptions axisopt = AxisOptions.create();
		axisopt.setIsLogScale(true);
		options.setHAxisOptions(axisopt);
		
		options.setTitleX(this.dataX);
		options.setTitleY(this.dataY);
		return options;
	}

	private AbstractDataTable createTable() {
		DataTable data = DataTable.create();
		data.removeRows(0, data.getNumberOfRows());
		data.addColumn(ColumnType.NUMBER, this.dataX);
		data.addColumn(ColumnType.NUMBER, this.dataY);
        for(Double d: this.values.keySet()){
        	int rowIndex = data.addRow();
        	data.setValue(rowIndex, 0, d);
        	data.setValue(rowIndex, 1, this.values.get(d));
        }		
		return data;
	}

	@Override
	public void setDefault() {
	}

}
