package amtrackservice.client.gui.elements;

import java.util.HashMap;

import amtrackservice.client.MapList;

import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.Widget;

public class AmTextField extends AmWidget  {
	
	private TextBox textbox;
	private String data;
	private MapList<String,String> preset;

	public AmTextField(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX) {
		super(label, datatype, description);
		textbox = new TextBox();
		this.preset = preset;
		this.data = dataX;
	}

	@Override
	public String getValue() {
		if (textbox.getText().equals("")){
			textbox.setText(preset.getValues().get(0));
		} return textbox.getValue();
	}

	@Override
	public Widget getWidget() {
		return textbox;
	}

	@Override
	public void setValue(HashMap<String, String> valueMap) {
		String str = valueMap.get(data);
		textbox.setValue(str);
	}
	
	@Override
	public String getDataLink() {
		return data;
	}


}
