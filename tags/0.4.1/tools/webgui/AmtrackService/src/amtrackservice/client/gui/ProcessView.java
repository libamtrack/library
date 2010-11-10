package amtrackservice.client.gui;

import amtrackservice.client.AmtrackService;

import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.ListBox;

public class ProcessView extends DialogBox {

	public ProcessView(final AmtrackService service, final String[] data) {
		// Set the dialog box's caption.
		setText("Dialog");

		// Enable animation.
		setAnimationEnabled(true);

		// Enable glass background.
		setGlassEnabled(true);
		
		final ListBox list = new ListBox(false);
		list.setVisibleItemCount(10);
		for(String s: data){
			list.addItem(s);
		}
		HorizontalPanel panel = new HorizontalPanel();

		panel.add(list);
		Button ok = new Button("OK");
		ok.addClickHandler(new ClickHandler() {
			public void onClick(ClickEvent event) {
//				service.openCalculation(Long.parseLong(data[list.getSelectedIndex()]));
				ProcessView.this.hide();
			}
		});
		panel.add(ok);
		setWidget(panel);
	}
}