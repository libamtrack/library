package amtrackservice.client.gui.elements;

import java.util.HashMap;

import amtrackservice.client.MapList;

import com.google.gwt.event.dom.client.KeyCodes;
import com.google.gwt.event.dom.client.KeyUpEvent;
import com.google.gwt.event.dom.client.KeyUpHandler;
import com.google.gwt.user.client.ui.DockPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.ListBox;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.Widget;

public class AmInputList extends AmWidget {

	private DockPanel widget;
	private ListBox listbox;
	private TextBox textbox;
	private String data;
	private MapList<String, String> preset;

	public AmInputList(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX) {
		super(label, datatype, description);
		widget = new DockPanel();
		listbox = new ListBox();
		textbox = new TextBox();
		this.preset = preset;
		
		for (String s : preset.getKeys()) {
			this.listbox.addItem(preset.getValue(s));
		}
		
		data = dataX;
		textbox.addKeyUpHandler(new KeyUpHandler() {

			@Override
			public void onKeyUp(KeyUpEvent event) {
				if (event.getNativeKeyCode() == KeyCodes.KEY_ENTER) {
					listbox.addItem(textbox.getText());
				}
			}
		});
		listbox.addKeyUpHandler(new KeyUpHandler() {

			@Override
			public void onKeyUp(KeyUpEvent event) {
				if (event.getNativeKeyCode() == KeyCodes.KEY_DELETE) {
					listbox.removeItem(listbox.getSelectedIndex());
				}
			}
		});
		listbox.setVisibleItemCount(5);

		listbox.setSize("100%", "100%");
		widget.add(listbox, DockPanel.CENTER);
		widget.add(textbox, DockPanel.SOUTH);
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
		return widget;
	}

	@Override
	public void setValue(HashMap<String, String> valueMap) {
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

}
