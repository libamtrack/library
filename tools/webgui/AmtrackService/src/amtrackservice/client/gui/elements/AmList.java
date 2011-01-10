package amtrackservice.client.gui.elements;

import java.util.HashMap;

import amtrackservice.client.MapList;

import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.Widget;

public class AmList extends AmWidget  {
	
	private ListBox listbox;
	private String data;
	private MapList<String,String> preset;

	public AmList(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX) {
		super(label, datatype, description);
		this.listbox = new ListBox();
		listbox.setVisibleItemCount(5);
		this.data = dataX;
		this.preset = preset;
	}

	@Override
	public String getValue() {
		String s = "";
		for (int i = 0; i < listbox.getItemCount(); i++) {
			s += " " + listbox.getItemText(i);
		}
		return s;
	}

	@Override
	public Widget getWidget() {
		// TODO Auto-generated method stub
		return listbox;
	}

	@Override
	public void setValue(HashMap<String, String> valueMap) {
		listbox.clear();
		String str = valueMap.get(data);
		if (str != null) {
			String[] s = str.split(" ");
			for (String c : s) {
				listbox.addItem(c);
			}
		}
		
	}
	
	@Override
	public String getDataLink() {
		return data;
	}

	@Override
	public void setDefault() {
		this.listbox.clear();
		for (String s : this.preset.getKeys()) {
			this.listbox.addItem(this.preset.getValue(s));
		}
	}

	@Override
	public void appendValue(HashMap<String, String> valueMap) {
		// TODO Auto-generated method stub
		
	}

}
