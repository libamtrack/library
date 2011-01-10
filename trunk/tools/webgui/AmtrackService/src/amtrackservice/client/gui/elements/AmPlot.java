package amtrackservice.client.gui.elements;

import java.util.HashMap;
import java.util.Random;
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
	private boolean xAxisLogByDefault;
	private boolean yAxisLogByDefault;
	
	private RadioButton xAxisScaleLogarithmicButton;		
	private RadioButton xAxisScaleLinearButton;

	private RadioButton yAxisScaleLogarithmicButton;		
	private RadioButton yAxisScaleLinearButton;	
	
	private HashMap< String, TreeMap<Double, Double>> values = new HashMap<String, TreeMap<Double,Double>>();

	public AmPlot(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX, String dataY,
			String dataZ, boolean xAxisLog, boolean yAxisLog) {
		super(label, datatype, description);

		Random rand = new Random(); 
		String id = label + rand.nextLong();

		this.dataX = dataX;
		this.dataY = dataY;
		this.xAxisLogByDefault = xAxisLog;
		this.yAxisLogByDefault = yAxisLog;
		
		xAxisScaleLogarithmicButton = new RadioButton("xaxis" + id, "logarithmic");		
		xAxisScaleLinearButton = new RadioButton("xaxis" + id, "linear");

		yAxisScaleLogarithmicButton = new RadioButton("yaxis" + id, "logarithmic");		
		yAxisScaleLinearButton = new RadioButton("yaxis" + id, "linear");	
		
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
		
		addDataSerie(valueMap, this.dataY);

		this.widget.clear();
		this.widget.add(createPlot());

        chart.draw(createTable(), createOptions());						
	}

	@Override
	public String getDataLink() {
		return null; // unimplemented
	}
			
	private Widget createPlot() {
        FlowPanel panel = new FlowPanel();        
        this.chart = new ScatterChart(DataTable.create(), createOptions() );
        		
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
        
        panel.clear();
		panel.add(axisScale);
        panel.add(this.chart);
        return panel;
	}		

	private Options createOptions() {
		Options options = Options.create();
		options.setWidth(640);
		options.setHeight(480);
		options.setTitle(this.getLabel().getText());		
		options.setLineSize(1);
		options.setPointSize(4);
		if( this.xAxisScaleLogarithmicButton.getValue() ){
			options.setTitleX("log( " + dataX + " )");
			this.xAxisScaleLogarithmicButton.setValue(true);
		} else {
			this.xAxisScaleLogarithmicButton.setValue(false);
			options.setTitleX(dataX);
		}
		
		if( this.yAxisScaleLogarithmicButton.getValue() ){
			this.yAxisScaleLogarithmicButton.setValue(true);
			options.setTitleY("log( " + dataY + " )");
		} else {
			this.yAxisScaleLogarithmicButton.setValue(false);
			options.setTitleY(dataY);
		}
		return options;
	}

	private AbstractDataTable createTable() {
		DataTable data = DataTable.create();
		data.removeRows(0, data.getNumberOfRows());
		data.addColumn(ColumnType.NUMBER, this.dataX);

		TreeMap<Double,Integer> xValuesIndexes = new TreeMap<Double, Integer>(); 
		for( String name : this.values.keySet()){
	        for(Double d: this.values.get(name).keySet()){
	        	xValuesIndexes.put(d, 0);
	        }
		}
		
		int i = 0;
		for( Double d: xValuesIndexes.keySet()){
			xValuesIndexes.put(d, i);
			i++;
		}
		
		data.addRows(xValuesIndexes.size());
		
		i = 1;
		for( String name : this.values.keySet()){
			data.addColumn(ColumnType.NUMBER, name);
	        for(Double d: this.values.get(name).keySet()){
	        	double xValueToInsert = 0., yValueToInsert = 0.;
	        	if( this.xAxisScaleLogarithmicButton.getValue() ){
	            	xValueToInsert = Math.log10(d);
	        	} else {
	        		xValueToInsert = d;
	        	}
	        	if( this.yAxisScaleLogarithmicButton.getValue() ){
	        		yValueToInsert = Math.log10(this.values.get(name).get(d));
	        	} else {
	        		yValueToInsert = this.values.get(name).get(d);
	        	}
	        	int rowIndex = xValuesIndexes.get(d);
            	data.setValue(rowIndex, 0, xValueToInsert);
        		data.setValue(rowIndex, i, yValueToInsert);
	        }
	        i++;
		}	
		
		return data;
	}

	public void setDefault() {
	}

	public void appendValue(HashMap<String, String> valueMap) {
		int newDataSerieNumber = this.values.size();

		addDataSerie(valueMap, this.dataY + "-" + newDataSerieNumber);
		
		this.widget.clear();		
		this.widget.add(createPlot());

        chart.draw(createTable(), createOptions());						
	}
	
	private void addDataSerie(HashMap<String, String> valueMap, String name){
		String[] xValues = {};
		if( valueMap.get(dataX) != null )
			xValues = valueMap.get(dataX).split(" ");
		String[] yValues = {};
		if( valueMap.get(dataY) != null )
			yValues = valueMap.get(dataY).split(" ");
		
		this.values.put(name, new TreeMap<Double, Double>());
		for (int i = 0; (i < xValues.length) && (i < yValues.length); i++) {
			double x = Double.parseDouble(xValues[i]);
			double y = Double.parseDouble(yValues[i]);
			this.values.get(name).put(x, y);
		}
	}

}
