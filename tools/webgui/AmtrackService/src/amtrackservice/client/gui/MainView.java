package amtrackservice.client.gui;

import amtrackservice.client.AmtrackService;
import amtrackservice.client.Logger;
import com.google.gwt.dom.client.Style.Unit;
import com.google.gwt.user.client.Command;
import com.google.gwt.user.client.ui.DockLayoutPanel;
import com.google.gwt.user.client.ui.MenuBar;
import com.google.gwt.user.client.ui.MenuItem;
import com.google.gwt.user.client.ui.RootLayoutPanel;
import com.google.gwt.user.client.ui.SplitLayoutPanel;
import com.google.gwt.user.client.ui.TabLayoutPanel;
import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.user.client.ui.TextArea;
import com.google.gwt.user.client.ui.Widget;

/**
 * 
 * This Represents the GUI 
 * @author Christoph Kolb, 2010 Hochschule Heilbronn
 */
public class MainView {

	private AmtrackService service;
	private TabLayoutPanel tabPanel;
	private TextArea console;
	private MenuBar menubar;
	private MenuBar menuNew;

	/**
	 * Constructor of MainView
	 * @param service the corresponding LibService
	 */
	public MainView(AmtrackService service) {
		this.service = service;
		Logger.init(this);
		console = new TextArea();
		
		menubar = new MenuBar();
		menuNew = new MenuBar(true);
		menubar.addItem("New calculation", menuNew);
		menubar.addSeparator();
		menubar.addItem("Close current calculation", new CalculationClose());

		tabPanel = new TabLayoutPanel(10, Unit.MM);
		
		SplitLayoutPanel rootPanel = new SplitLayoutPanel();
		//rootPanel.addSouth(console,200);
		//console.setSize("100%", "90%");
		DockLayoutPanel main = new DockLayoutPanel(Unit.MM);
		main.addNorth(menubar, 5);
		main.add(tabPanel);
		
		rootPanel.add(main);
		
		RootLayoutPanel rp  = RootLayoutPanel.get();
		rp.add(rootPanel);
		
	}
	
	/**
	 * Adds a new Widget to the Central Panel.
	 * @param input The widget
	 * @param name A title for the widget
	 */
	public void addTabPanel(Widget input, String name){
		tabPanel.add(input,name);
		tabPanel.selectTab(input);		
	}
	
	/**
	 * 
	 * @param line
	 */
	public void addConsoleLine(String line){
		console.setText(console.getText()+"\n"+line);
		console.getElement().setScrollTop(console.getElement().getScrollHeight());
	}
	
	private class CalculationStart implements Command{

		private String method;
		
		public CalculationStart(String method) {
			this.method = method;
		}
		@Override
		public void execute() {
			service.openCalculation(method);
			
		}
		
	}
	
	private class CalculationClose implements Command{
		
		@Override
		public void execute() {
			int index = tabPanel.getSelectedIndex();
			tabPanel.remove(index);
		}
		
	}

	public void initMenu(String[] array) {
		for(String m : array){
			MenuItem mi = new MenuItem(m, new CalculationStart(m));
			menuNew.addItem(mi);
		}
		
		
	}
}