package amtrackservice.client.gui.elements;

import java.util.HashMap;

import amtrackservice.client.MapList;
import amtrackservice.client.gui.elements.ScatterChartLogAxis.AxisOptions;
import amtrackservice.client.gui.elements.ScatterChartLogAxis.Options;

import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.ui.FlowPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.ToggleButton;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.AbstractDataTable.ColumnType;
import com.google.gwt.visualization.client.DataTable;

public class AmPlot extends AmWidget {

	private String dataX;
	private String dataY;
	private FlowPanel widget = new FlowPanel();
	private ScatterChartLogAxis chart;
	
	private HashMap<Double, Double> values = new HashMap<Double, Double>();

	public AmPlot(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX, String dataY,
			String dataZ) {
		super(label, datatype, description);
		this.dataX = dataX;
		this.dataY = dataY;
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
		final ToggleButton tbxaxis = new ToggleButton("switch X axis to logarithmic scale","switch X axis to linear scale");
		tbxaxis.addClickHandler(new ClickHandler() {			
			public void onClick(ClickEvent event) {
				if( tbxaxis.isDown() ){
					setXAxisLogScale(true);
				} else {
					setXAxisLogScale(false);
				}
			}
		});
		widget.add(tbxaxis);
		final ToggleButton tbyaxis = new ToggleButton("switch Y axis to logarithmic scale","switch Y axis to linear scale");
		tbyaxis.addClickHandler(new ClickHandler() {			
			public void onClick(ClickEvent event) {
				if( tbyaxis.isDown() ){
					setYAxisLogScale(true);
				} else {
					setYAxisLogScale(false);
				}
			}
		});
		widget.add(tbyaxis);
		this.widget.add(createPlot());
	}

	@Override
	public String getDataLink() {
		return null; // unimplemented
	}
	
	public void setXAxisLogScale(boolean logscale){
		Options options = createOptions();
		AxisOptions axisopt = AxisOptions.create();
		axisopt.setIsLogScale(logscale);
		options.setHAxisOptions(axisopt);
		this.chart.draw(createTable(), options);
	}
	
	public void setYAxisLogScale(boolean logscale){
		Options options = createOptions();
		AxisOptions axisopt = AxisOptions.create();
		axisopt.setIsLogScale(logscale);
		options.setVAxisOptions(axisopt);
		this.chart.draw(createTable(), options);
	}
		
	private Widget createPlot() {
        FlowPanel panel = new FlowPanel();        
        this.chart = new ScatterChartLogAxis(createTable(), createOptions() );
        panel.clear();
        panel.add(this.chart);
        return panel;

	}		

	private Options createOptions() {
		Options options = Options.create();
		options.setWidth(640);
		options.setHeight(480);
		options.setTitle(this.getLabel().getText());
		options.setLineSize(1);
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
