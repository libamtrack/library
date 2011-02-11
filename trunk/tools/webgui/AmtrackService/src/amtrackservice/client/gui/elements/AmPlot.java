package amtrackservice.client.gui.elements;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
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


/**
 *  one data serie, which contains:
 *    vector of X and corresponding to them Y values
 *    names (types) of X and Y data
 *    name of serie 
 **/
class DataSerie {

	private String name;
	private String dataXname;
	private String dataYname;
	private TreeMap<Double, Double> dataXY;
	
	/**
	 * @param name
	 * @param dataXY
	 */
	public DataSerie(String name, String dataXname, String dataYname, TreeMap<Double, Double> dataXY) {
		this.name = name;
		this.dataXname = dataXname;
		this.dataYname = dataYname;
		this.dataXY = dataXY;
	}		
	
	/**
	 * TODO
	 * @return
	 */
	public int size(){
		return dataXY.size();
	}
	
	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}
	
	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * @param dataXY the dataXY to set
	 */
	public void setData(TreeMap<Double, Double> dataXY) {
		this.dataXY = dataXY;
	}
	
	/**
	 * @return the dataXY
	 */
	public double get(double x) {
		return dataXY.get(x);
	}
	
	/**
	 * TODO
	 * @return
	 */
	public Collection<Double> getDataX(){
		return dataXY.keySet();
	}
	
	/**
	 * TODO
	 * @return
	 */
	public Collection<Double> getDataY(){
		return dataXY.values();
	}

	/**
	 * @param dataXname the dataXname to set
	 */
	public void setDataXname(String dataXname) {
		this.dataXname = dataXname;
	}

	/**
	 * @return the dataXname
	 */
	public String getDataXname() {
		return dataXname;
	}

	/**
	 * @param dataYname the dataYname to set
	 */
	public void setDataYname(String dataYname) {
		this.dataYname = dataYname;
	}

	/**
	 * @return the dataYname
	 */
	public String getDataYname() {
		return dataYname;
	}
		
	/**
	 * TODO
	 * @param xAxisLogScale
	 * @param yAxisLogScale
	 * @return
	 */
	public AbstractDataTable getDataTable( boolean xAxisLogScale , boolean yAxisLogScale){
		DataTable data = DataTable.create();
		data.removeRows(0, data.getNumberOfRows());
		data.addColumn(ColumnType.NUMBER, this.getDataXname());
		data.addColumn(ColumnType.NUMBER, this.getDataYname());
		data.addRows(this.dataXY.size());
				
		int rowIndex = 0;
        for( Map.Entry<Double, Double> item : this.dataXY.entrySet()){        
        	double xValueToInsert = 0., yValueToInsert = 0.;
        	if( xAxisLogScale ){
            	xValueToInsert = Math.log10(item.getKey());
        	} else {
        		xValueToInsert = item.getKey();
        	}
        	if( yAxisLogScale ){
        		yValueToInsert = Math.log10(item.getValue());
        	} else {
        		yValueToInsert = item.getValue();
        	}        	
           	data.setValue(rowIndex, 0, xValueToInsert);
       		data.setValue(rowIndex, 1, yValueToInsert);
       		rowIndex++;
        }
		
		return data;
	}
}

/**
 *  collection of data series with the same name (type) of data
 *  at X and Y axis, contains:
 *    collection of DataSeries
 *    names of X and Y data
 *    name of collection 
 **/
class DataSerieCollection {

	private String name;
	private String dataXname;
	private String dataYname;
	private TreeMap<String, DataSerie> dataSerieMap;
	
	/**
	 * TODO
	 * @param name
	 */
    public DataSerieCollection(String name) {
        this.setName(name);
        dataSerieMap = new TreeMap<String, DataSerie>();
    }

    /**
     * TODO
     * @return
     */
    public int size(){
    	return dataSerieMap.size();
    }
    
    /**
     * TODO
     */
    public void clear(){
    	this.dataSerieMap.clear();
    }
    
	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the dataXname
	 */
	public String getDataXname() {
		return dataXname;
	}

	/**
	 * @return the dataYname
	 */
	public String getDataYname() {
		return dataYname;
	}
	
	/**
	 * TODO
	 * @param serie
	 * @return
	 */
	public int add(DataSerie serie){
		if( this.size() == 0 ){
			this.dataSerieMap.put(serie.getName(), serie);
			this.dataXname = serie.getDataXname();
			this.dataYname = serie.getDataYname();
		} else {
			/** problem - data X name in collection is different from data X name in serie */
			if( !this.getDataXname().equals(serie.getDataXname()) ){
				return -1;
			}
			/** problem - data Y name in collection is different from data Y name in serie */
			if( !this.getDataYname().equals(serie.getDataYname()) ){
				return -2;
			}
			/** problem - another data serie with the same name exists in collection */
			if( this.dataSerieMap.containsKey(serie.getName())){
				return -3;
			}
			this.dataSerieMap.put(serie.getName(), serie);
		}
		return 0;
	}
	
	
	/**
	 * TODO
	 * @return
	 */
	private boolean checkIfAllSeriesHaveTheSameDataX(){
		if( !this.dataSerieMap.isEmpty()){
			// first dataX item:
			String firstKey = this.dataSerieMap.firstKey();
			Collection<Double> xVector = this.dataSerieMap.get(firstKey).getDataX();
			// loop over all items
			for( Map.Entry<String, DataSerie> item : this.dataSerieMap.entrySet()){
				// if found one item which is different than the first one - exit with false
				if( !item.getValue().getDataX().equals(xVector))
					return false;
			}
		}
		return true;
	}
	
	
	/**
	 * TODO
	 * @param flag
	 * @param value
	 * @return
	 */
	private double convertToLogarithm(boolean flag, double value){
		if( flag )
			return Math.log10(value);
		return value;		
	}
	
	
	/**
	 * TODO
	 * @param xAxisLogScale
	 * @param yAxisLogScale
	 * @return
	 */
	public AbstractDataTable getDataTable( boolean xAxisLogScale , boolean yAxisLogScale){
		DataTable data = DataTable.create();
		
		// clean datatable
		data.removeRows(0, data.getNumberOfRows());
		data.removeColumns(0, data.getNumberOfColumns());
		
		// add first column for storing X values 
		data.addColumn(ColumnType.NUMBER, this.getDataXname());
		
		boolean allDataXVectorsTheSame = checkIfAllSeriesHaveTheSameDataX();
		
		if( allDataXVectorsTheSame ){
			// get first data X from first data serie 		
			String firstKey = this.dataSerieMap.firstKey();
			Collection<Double> xVector = this.dataSerieMap.get(firstKey).getDataX();
			
			// add sufficient amount of rows
			data.addRows(xVector.size());
			
			// add sufficient amount of columns
			for( Map.Entry<String, DataSerie> item : this.dataSerieMap.entrySet()){
				data.addColumn(ColumnType.NUMBER, item.getKey());
			}
			
			// loop over all x from data X
			int rowIndex = 0;
			for( Double x : xVector){
				
				// insert X into a data table
	        	double xValueToInsert = convertToLogarithm(xAxisLogScale, x);
				data.setValue(rowIndex, 0, xValueToInsert);
				
				// loop over all data series and insert Y value corresponding to inserted X
				int columnIndex = 1;
				for( Map.Entry<String, DataSerie> item : this.dataSerieMap.entrySet()){
		        	double yValueToInsert = convertToLogarithm(yAxisLogScale, item.getValue().get(x));
					data.setValue(rowIndex, columnIndex, yValueToInsert);
					columnIndex++;
				}
				rowIndex++;
			}
		} else {
			int columnIndex = 1;
			int rowIndex = 0;
			// loop over all series
			for( Map.Entry<String, DataSerie> item : this.dataSerieMap.entrySet()){
				// for every serie append some amount of rows and columns
				data.addRows(item.getValue().size());
				data.addColumn(ColumnType.NUMBER, item.getKey());
				
				// loop over all x from data X
				for( Double x : item.getValue().getDataX()){
		        	
					// insert X
					double xValueToInsert = convertToLogarithm(xAxisLogScale, x);
					data.setValue(rowIndex, 0, xValueToInsert);
					
					// insert Y
		        	double yValueToInsert = convertToLogarithm(yAxisLogScale, item.getValue().get(x));			
					data.setValue(rowIndex, columnIndex, yValueToInsert);						
					
					rowIndex++;					
				}
				columnIndex++;
			}
		}
		return data;
	}
	
}


/**
 * TODO
 * @author grzanka
 */
public class AmPlot extends AmWidget {
	
	private String dataX;
	private String dataY;
	private String labelX;
	private String labelY;
	private FlowPanel widget = new FlowPanel();
	private ScatterChart chart;
	private boolean xAxisLogByDefault;
	private boolean yAxisLogByDefault;
	
	private RadioButton xAxisScaleLogarithmicButton;		
	private RadioButton xAxisScaleLinearButton;

	private RadioButton yAxisScaleLogarithmicButton;		
	private RadioButton yAxisScaleLinearButton;	
	
	private DataSerieCollection dataSerieCollection;

	/**
	 * TODO
	 * @param plotLabel
	 * @param datatype
	 * @param description
	 * @param preset
	 * @param dataX
	 * @param dataY
	 * @param labelX
	 * @param labelY
	 * @param xAxisLog
	 * @param yAxisLog
	 */
	public AmPlot(String plotLabel, String datatype, HTML description,
			MapList<String, String> preset, 
			String dataX, String dataY,
			String labelX, String labelY,
			boolean xAxisLog, boolean yAxisLog) {
		super(plotLabel, datatype, description);

		Random rand = new Random(); 
		String id = plotLabel + rand.nextLong();

		this.dataX = dataX;
		this.dataY = dataY;
		this.labelX = labelX;
		this.labelY = labelY;
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
		
		this.dataSerieCollection = new DataSerieCollection(plotLabel);
		
	}

	
	/**
	 * TODO
	 */
	public Widget getWidget() {
		return widget;
	}

	
	/**
	 * TODO
	 */
	public int setValue(HashMap<String, String> valueMap) {		
		this.dataSerieCollection.clear();		
		int status = this.appendValue(valueMap);
		return status;
	}
			
	
	/**
	 * TODO
	 * @return
	 */
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

	
	/**
	 * TODO
	 * @return
	 */
	private Options createOptions() {
		Options options = Options.create();
		options.setWidth(640);
		options.setHeight(480);
		options.setTitle(this.getLabel().getText());		
		options.setLineSize(1);
		options.setPointSize(4);
		if( this.xAxisScaleLogarithmicButton.getValue() ){
			options.setTitleX("log( " + labelX + " )");
			this.xAxisScaleLogarithmicButton.setValue(true);
		} else {
			this.xAxisScaleLogarithmicButton.setValue(false);
			options.setTitleX(labelX);
		}
		
		if( this.yAxisScaleLogarithmicButton.getValue() ){
			this.yAxisScaleLogarithmicButton.setValue(true);
			options.setTitleY("log( " + labelY + " )");
		} else {
			this.yAxisScaleLogarithmicButton.setValue(false);
			options.setTitleY(labelY);
		}
		return options;
	}

	/**
	 * TODO
	 * @return
	 */
	private AbstractDataTable createTable() {
		return this.dataSerieCollection.getDataTable(this.xAxisScaleLogarithmicButton.getValue(), yAxisScaleLogarithmicButton.getValue());
	}

	/**
	 * TODO
	 */
	public int appendValue(HashMap<String, String> valueMap) {
		
		String serieLabel = valueMap.get("label");
				
		String[] xValues = {};
		if( valueMap.get(dataX) != null )
			xValues = valueMap.get(dataX).split(" ");
		String[] yValues = {};
		if( valueMap.get(dataY) != null )
			yValues = valueMap.get(dataY).split(" ");
		
		TreeMap<Double, Double> dataXY = new TreeMap<Double, Double>();
		for (int i = 0; (i < xValues.length) && (i < yValues.length); i++) {
			double x = Double.parseDouble(xValues[i]);
			double y = Double.parseDouble(yValues[i]);
			dataXY.put(x, y);
		}
				
		DataSerie serie = new DataSerie(serieLabel, this.dataX, this.dataY, dataXY);		
		int status = this.dataSerieCollection.add(serie);
		
		if( status == 0){		
			this.widget.clear();
			this.widget.add(createPlot());
			AbstractDataTable dataTable = createTable();
			if( dataTable == null ){
				return -4;
			}
			chart.draw(dataTable, createOptions());
		}
		return status;
	}
	
	/**
	 * TODO
	 */
	public String getValue() {
		return null;
	}

	/**
	 * TODO
	 */
	public void setDefault() {}

	/**
	 * TODO
	 */
	public String getDataLink() {
		return null;
	}

}
