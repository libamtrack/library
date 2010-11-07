package amtrackservice.client.gui.elements;

import java.util.HashMap;
import java.util.TreeMap;

import amtrackservice.client.MapList;
import amtrackservice.client.gui.elements.ScatterChartLogAxis.AxisOptions;
import amtrackservice.client.gui.elements.ScatterChartLogAxis.Options;

import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.ui.FlowPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.HasVerticalAlignment;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.RadioButton;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.AbstractDataTable.ColumnType;
import com.google.gwt.visualization.client.DataTable;

public class AmPlot extends AmWidget {

	private String dataX;
	private String dataY;
	private FlowPanel widget = new FlowPanel();
	private ScatterChartLogAxis chart;
	boolean xAxisLog;
	boolean yAxisLog;
	
	private TreeMap<Double, Double> values = new TreeMap<Double, Double>();

	public AmPlot(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX, String dataY,
			String dataZ, boolean xAxisLog, boolean yAxisLog) {
		super(label, datatype, description);
		this.dataX = dataX;
		this.dataY = dataY;
		this.xAxisLog = xAxisLog;
		this.yAxisLog = yAxisLog;
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
				
		VerticalPanel axisScale = new VerticalPanel();
		axisScale.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_LEFT);
		
		HorizontalPanel xAxisScale = new HorizontalPanel();
		xAxisScale.setVerticalAlignment(HasVerticalAlignment.ALIGN_MIDDLE);

		HorizontalPanel yAxisScale = new HorizontalPanel();
		yAxisScale.setVerticalAlignment(HasVerticalAlignment.ALIGN_MIDDLE);

		
		HTML xAxisButtonLabel = new HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;X axis scale: ");
		
		final RadioButton xAxisScaleLogarithmicButton = new RadioButton("xaxis", "logarithmic");		
		final RadioButton xAxisScaleLinearButton = new RadioButton("xaxis", "linear");
		if( xAxisLog ){
			xAxisScaleLogarithmicButton.setValue(true);
		} else {
			xAxisScaleLinearButton.setValue(true);			
		}
		HTML yAxisButtonLabel = new HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Y axis scale: ");
		
		final RadioButton yAxisScaleLogarithmicButton = new RadioButton("yaxis", "logarithmic");		
		final RadioButton yAxisScaleLinearButton = new RadioButton("yaxis", "linear");	
		if( yAxisLog ){
			yAxisScaleLogarithmicButton.setValue(true);
		} else {
			yAxisScaleLinearButton.setValue(true);			
		}
		
		xAxisScale.add(xAxisButtonLabel);
		xAxisScale.add(xAxisScaleLinearButton);
		xAxisScale.add(xAxisScaleLogarithmicButton);
		yAxisScale.add(yAxisButtonLabel);
		yAxisScale.add(yAxisScaleLinearButton);
		yAxisScale.add(yAxisScaleLogarithmicButton);
		
		axisScale.add(xAxisScale);
		axisScale.add(yAxisScale);
		
		ClickHandler xAxisLinClickHandler = new ClickHandler() {			
			public void onClick(ClickEvent event) {
				setAxisLogScales(false, yAxisScaleLogarithmicButton.getValue());
			}
		};
		ClickHandler xAxisLogClickHandler = new ClickHandler() {			
			public void onClick(ClickEvent event) {
				setAxisLogScales(true, yAxisScaleLogarithmicButton.getValue());
			}
		};
		ClickHandler yAxisLinClickHandler = new ClickHandler() {			
			public void onClick(ClickEvent event) {
				setAxisLogScales(xAxisScaleLogarithmicButton.getValue(), false);
			}
		};
		ClickHandler yAxisLogClickHandler = new ClickHandler() {			
			public void onClick(ClickEvent event) {
				setAxisLogScales(xAxisScaleLogarithmicButton.getValue(), true);
			}
		};
		xAxisScaleLinearButton.addClickHandler(xAxisLinClickHandler);
		xAxisScaleLogarithmicButton.addClickHandler(xAxisLogClickHandler);
		yAxisScaleLinearButton.addClickHandler(yAxisLinClickHandler);
		yAxisScaleLogarithmicButton.addClickHandler(yAxisLogClickHandler);
		
		widget.add(axisScale);
		this.widget.add(createPlot());
		setAxisLogScales(xAxisScaleLogarithmicButton.getValue(), yAxisScaleLogarithmicButton.getValue());
	}

	@Override
	public String getDataLink() {
		return null; // unimplemented
	}
	
	public void setAxisLogScales(boolean logscaleX, boolean logscaleY){
		Options options = createOptions();
		AxisOptions axisoptX = AxisOptions.create();
		AxisOptions axisoptY = AxisOptions.create();
		axisoptX.setIsLogScale(logscaleX);
		options.setHAxisOptions(axisoptX);
		axisoptY.setIsLogScale(logscaleY);
		options.setVAxisOptions(axisoptY);
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
