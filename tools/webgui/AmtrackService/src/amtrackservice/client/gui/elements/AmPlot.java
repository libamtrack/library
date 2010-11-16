package amtrackservice.client.gui.elements;

import java.util.HashMap;
import java.util.TreeMap;

import amtrackservice.client.MapList;

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
import com.google.gwt.visualization.client.visualizations.ScatterChart;
import com.google.gwt.visualization.client.visualizations.ScatterChart.Options;

public class AmPlot extends AmWidget {

	private String dataX;
	private String dataY;
	private FlowPanel widget = new FlowPanel();
	private ScatterChart chart;
	boolean xAxisLogByDefault;
	boolean yAxisLogByDefault;
	
	RadioButton xAxisScaleLogarithmicButton = new RadioButton("xaxis", "logarithmic");		
	RadioButton xAxisScaleLinearButton = new RadioButton("xaxis", "linear");

	RadioButton yAxisScaleLogarithmicButton = new RadioButton("yaxis", "logarithmic");		
	RadioButton yAxisScaleLinearButton = new RadioButton("yaxis", "linear");	
	
	private TreeMap<Double, Double> values = new TreeMap<Double, Double>();

	public AmPlot(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX, String dataY,
			String dataZ, boolean xAxisLog, boolean yAxisLog) {
		super(label, datatype, description);
		this.dataX = dataX;
		this.dataY = dataY;
		this.xAxisLogByDefault = xAxisLog;
		this.yAxisLogByDefault = yAxisLog;
		
		xAxisScaleLogarithmicButton.setValue(xAxisLogByDefault);
		yAxisScaleLogarithmicButton.setValue(yAxisLogByDefault);
		xAxisScaleLinearButton.setValue(!xAxisLogByDefault);
		yAxisScaleLinearButton.setValue(!yAxisLogByDefault);
		
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
		HTML yAxisButtonLabel = new HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Y axis scale: ");
				
		xAxisScale.add(xAxisButtonLabel);
		xAxisScale.add(this.xAxisScaleLinearButton);
		xAxisScale.add(this.xAxisScaleLogarithmicButton);
		yAxisScale.add(yAxisButtonLabel);
		yAxisScale.add(this.yAxisScaleLinearButton);
		yAxisScale.add(this.yAxisScaleLogarithmicButton);
		
		axisScale.add(xAxisScale);
		axisScale.add(yAxisScale);
		
		// one click handler for all buttons which checks state 

		ClickHandler generalAxisHandler = new ClickHandler() {			
			public void onClick(ClickEvent event) {
				chart.draw(createTable(), createOptions());
			}
		};

		xAxisScaleLinearButton.addClickHandler(generalAxisHandler);
		xAxisScaleLogarithmicButton.addClickHandler(generalAxisHandler);
		yAxisScaleLinearButton.addClickHandler(generalAxisHandler);
		yAxisScaleLogarithmicButton.addClickHandler(generalAxisHandler);

		widget.add(axisScale);
		this.widget.add(createPlot());
		
        chart.draw(createTable(), createOptions());						
	}

	@Override
	public String getDataLink() {
		return null; // unimplemented
	}
			
	private Widget createPlot() {
        FlowPanel panel = new FlowPanel();        
        this.chart = new ScatterChart(createTable(), createOptions() );
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
		if( xAxisScaleLogarithmicButton.getValue() ){
			options.setTitleX("log( " + dataX + " )");
		} else {
			options.setTitleX(dataX);
		}
		
		if( yAxisScaleLogarithmicButton.getValue() ){
			options.setTitleY("log( " + dataY + " )");
		} else {
			options.setTitleY(dataY);
		}
		return options;
	}

	private AbstractDataTable createTable() {
		DataTable data = DataTable.create();
		data.removeRows(0, data.getNumberOfRows());
		data.addColumn(ColumnType.NUMBER, this.dataX);
		data.addColumn(ColumnType.NUMBER, this.dataY);		
        for(Double d: this.values.keySet()){
        	int rowIndex = data.addRow();
        	if( this.xAxisScaleLogarithmicButton.getValue() ){
            	data.setValue(rowIndex, 0, Math.log10(d));        		
        	} else {
            	data.setValue(rowIndex, 0, d);
        	}
        	if( this.yAxisScaleLogarithmicButton.getValue() ){
        		data.setValue(rowIndex, 1, Math.log10(this.values.get(d)));
        	} else {
        		data.setValue(rowIndex, 1, this.values.get(d));
        	}
        }		
		return data;
	}

	@Override
	public void setDefault() {
	}

}
