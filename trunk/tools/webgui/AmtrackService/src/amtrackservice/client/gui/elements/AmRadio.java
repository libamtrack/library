package amtrackservice.client.gui.elements;

import java.util.HashMap;
import amtrackservice.client.MapList;

import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.RadioButton;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;

public class AmRadio extends AmWidget {

	private VerticalPanel panelWithRadioButtons;
	private MapList<String,String> preset;
	private String data;
	MapList<String, String> entry;
	
	public AmRadio(String label, String datatype, HTML description,
			MapList<String, String> entry, MapList<String, String> preset,
			String dataX) {
		super(label, datatype, description);
		this.preset = preset;
		this.entry = entry;
		this.panelWithRadioButtons = new VerticalPanel();
		
		// create group of radio buttons under main widget name (label)
		for(String item : entry.getKeys()){
			String value = entry.getValue(item);
			RadioButton button = new RadioButton(label,item);
			button.setValue(false);
			if(value.equals(preset.getValues().get(0)))
				button.setValue(true);
			this.panelWithRadioButtons.add(button);
		}
		this.data = dataX;
	}

	@Override
	public String getValue() {
		String result = "";
		for( int index = 0; index < this.panelWithRadioButtons.getWidgetCount() ; index++){
			RadioButton button = (RadioButton) this.panelWithRadioButtons.getWidget(index); 
			if( button.getValue() ){
				String value = button.getText();
				result = this.entry.getValue(value);
			}
		}
		return result;
	}

	@Override
	public Widget getWidget() {
		return panelWithRadioButtons;
	}

	@Override
	public void setValue(HashMap<String, String> valueMap) {
		// TODO
	}

	@Override
	public String getDataLink() {
		return data;
	}

	@Override
	public void setDefault() {
		for( int index = 0; index < this.panelWithRadioButtons.getWidgetCount() ; index++){
			RadioButton button = (RadioButton) this.panelWithRadioButtons.getWidget(index);
			button.setValue(false);
			String name = button.getText(); 
			String value = this.entry.getValue(name);
			if( value.equals(preset.getValues().get(0)) ){
				button.setValue(true);
			}
		}			
	}

	@Override
	public void appendValue(HashMap<String, String> valueMap) {
		// TODO Auto-generated method stub
		
	}

}
