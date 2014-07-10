package amtrackservice.client.gui.elements;

import java.util.HashMap;

import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.Label;
import com.google.gwt.user.client.ui.Widget;

public abstract class AmWidget {
	private HTML description;
	private String datatype;
	private String label;
	private boolean canHaveMultipleDataSeries;

	public AmWidget(String label, String datatype, HTML description) {
		this.description = description;
		this.label = label;
		this.datatype = datatype;
		this.canHaveMultipleDataSeries = false;
	}
	public abstract String getValue();

	public abstract int setValue(HashMap<String, String> valueMap);
	
	public abstract int appendValue(HashMap<String, String> valueMap);
	
	public abstract void setDefault();

	public HTML getDescription() {
		return description;
	}

	public void setDescription(HTML description) {
		this.description = description;
	}

	public String getDatatype() {
		return datatype;
	}
	
	public abstract String getDataLink();

	public Label getLabel() {
		return new Label(label);
	}

	public abstract Widget getWidget();
	
	public void setMultipleDataSeriesEnable( boolean flag ){
		this.canHaveMultipleDataSeries = flag;
	}
	
	public boolean getMultipleDataSeriesEnable(){
		return this.canHaveMultipleDataSeries;
	}
	
}
