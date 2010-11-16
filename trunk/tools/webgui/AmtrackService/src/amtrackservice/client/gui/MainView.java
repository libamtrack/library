package amtrackservice.client.gui;

import amtrackservice.client.AmtrackService;
import amtrackservice.client.AmtrackService.AmtrackServiceResources;
import amtrackservice.client.Logger;

import com.google.gwt.core.client.GWT;
import com.google.gwt.dom.client.Style.Unit;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.ui.DecoratedStackPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.HasVerticalAlignment;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.PushButton;
import com.google.gwt.user.client.ui.RootLayoutPanel;
import com.google.gwt.user.client.ui.TabLayoutPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;

/**
 * This Represents the GUI 
 * @author Christoph Kolb, 2010 Hochschule Heilbronn
 */
public class MainView {

	private TabLayoutPanel centralTabPanel;
	private DecoratedStackPanel leftDownPanelFunctions;
	private HorizontalPanel upperRightPanelCloseCurrentTab;

	/**
	 * Constructor of MainView
	 * @param service 
	 */
	public MainView(AmtrackService service) {
		Logger.init(this);
		
		AmtrackServiceResources resources = GWT.create(AmtrackServiceResources.class);
		
		/**  Constructors **/
		RootLayoutPanel rootPanel  = RootLayoutPanel.get();
	    HorizontalPanel wholePagePanel = new HorizontalPanel();
		VerticalPanel leftPanel = new VerticalPanel();
	    VerticalPanel rightPanel = new VerticalPanel();

	    Image leftPanelUpperLogo = new Image(resources.logo());
		leftDownPanelFunctions = new DecoratedStackPanel();
	    
		final HTML rightMainTab = new HTML("<center><big><font size=\"+2\"><br><br>Welcome to libamtrack webGUI</font></big><br><br>" + 
				"<big>HOWTO</big><br><br>" +
				"<font size=\"1\">1. Choose a function below.<br>"+
				"2. Enter parameters and recalculate.<br>" +
				"3. Results printed on the right.<br>"+
				"4. Open another/close function.<br>"+
                "Please report any problems/feedback<br>"+
                "to: <b><u><a href=\"mailto:Leszek.Grzanka@ifj.edu.pl\" target=\"_blank\">Leszek.Grzanka@ifj.edu.pl</a></u></b><br></font></center>");
		
		upperRightPanelCloseCurrentTab = new HorizontalPanel();
	    HTML htmlCloseCurrentTab = new HTML("Close current calculation: ", true);
	    Image imageCloseCurrentTab = new Image(resources.close());
	    PushButton closeCurrentTabButton = new PushButton(imageCloseCurrentTab);
		
	    centralTabPanel = new TabLayoutPanel(10, Unit.MM);

		/**  Tree of life **/
		rootPanel.add(wholePagePanel);		
		wholePagePanel.add(leftPanel);		
	    wholePagePanel.add(rightPanel);	    
		
	    leftPanel.add(leftPanelUpperLogo);
	    leftPanel.add(leftDownPanelFunctions);

	    rightPanel.add(upperRightPanelCloseCurrentTab);
	    upperRightPanelCloseCurrentTab.add(htmlCloseCurrentTab);
	    upperRightPanelCloseCurrentTab.add(closeCurrentTabButton);
	    rightPanel.add(centralTabPanel);
	    
		
		/**  Dimensions **/
	    rootPanel.setSize("1180px", "800px");
		rootPanel.setWidgetLeftWidth(wholePagePanel, 0.0, Unit.PX, 1180.0, Unit.PX);
		rootPanel.setWidgetTopHeight(wholePagePanel, 0.0, Unit.PX, 800.0, Unit.PX);

		leftPanel.setWidth("180px");
	    leftPanel.setCellHeight(leftPanelUpperLogo, "120px");
	    leftPanel.setCellWidth(leftDownPanelFunctions, "180px");
		
	    rightPanel.setWidth("600px");
	    rightPanel.setCellWidth(upperRightPanelCloseCurrentTab, "600px");

	    leftPanelUpperLogo.setSize("180px", "120px");
	    leftDownPanelFunctions.setWidth("180px");
	    
	    upperRightPanelCloseCurrentTab.setSize("400px", "32px");
	    upperRightPanelCloseCurrentTab.setCellHeight(htmlCloseCurrentTab, "32px");
	    upperRightPanelCloseCurrentTab.setCellHeight(closeCurrentTabButton, "32px");
	    upperRightPanelCloseCurrentTab.setCellWidth(closeCurrentTabButton, "32px");
	    htmlCloseCurrentTab.setSize("200px", "32px");
	    closeCurrentTabButton.setSize("32px", "32px");
	    imageCloseCurrentTab.setSize("32px", "32px");

	    centralTabPanel.setSize("1000px", "800px");
		
		
		/**  Alignment **/
		wholePagePanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_CENTER);

		leftPanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_CENTER);
	    leftPanel.setCellVerticalAlignment(leftPanelUpperLogo, HasVerticalAlignment.ALIGN_BOTTOM);
	    leftPanel.setCellHorizontalAlignment(leftPanelUpperLogo, HasHorizontalAlignment.ALIGN_CENTER);
	    leftPanel.setCellHorizontalAlignment(leftDownPanelFunctions, HasHorizontalAlignment.ALIGN_CENTER);
	    
	    rightPanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_LEFT);
	    rightPanel.setCellHorizontalAlignment(upperRightPanelCloseCurrentTab, HasHorizontalAlignment.ALIGN_LEFT);
	    upperRightPanelCloseCurrentTab.setCellHorizontalAlignment(htmlCloseCurrentTab, HasHorizontalAlignment.ALIGN_RIGHT);
	    upperRightPanelCloseCurrentTab.setCellVerticalAlignment(htmlCloseCurrentTab, HasVerticalAlignment.ALIGN_MIDDLE);
	    upperRightPanelCloseCurrentTab.setCellHorizontalAlignment(closeCurrentTabButton, HasHorizontalAlignment.ALIGN_RIGHT);
	    upperRightPanelCloseCurrentTab.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_RIGHT);
	    htmlCloseCurrentTab.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_RIGHT);		
		
	    
		/**  Others **/
	    ClickHandler closeHandler = new ClickHandler() {					
			@Override
			public void onClick(ClickEvent event) {
				int index = centralTabPanel.getSelectedIndex();
				if( index > -1 ){
					centralTabPanel.remove(index);
				}
				if( centralTabPanel.getWidgetCount() == 0){
					upperRightPanelCloseCurrentTab.setVisible(false);
					centralTabPanel.add(rightMainTab);
				}
			}
		};
	    closeCurrentTabButton.addClickHandler(closeHandler);
	    htmlCloseCurrentTab.setStyleName("h1");
	    upperRightPanelCloseCurrentTab.setVisible(false);

		if( centralTabPanel.getWidgetCount() == 0){
			centralTabPanel.add(rightMainTab);
		}

	}
	
	/**
	 * Adds a new Widget to the Central Panel.
	 * @param input The widget
	 * @param name A title for the widget
	 */
	public void addTabPanel(Widget input, String name){
		if( centralTabPanel.getWidgetCount() == 1){
			Widget w = centralTabPanel.getWidget(0);
			if( w.getClass().getName().equals("com.google.gwt.user.client.ui.HTML") ){
				centralTabPanel.remove(0);
			}
		}		

		centralTabPanel.add(input,name);
		centralTabPanel.selectTab(input);
		upperRightPanelCloseCurrentTab.setVisible(true);
	}
	
	public void addWidgetToLeftPanel(Widget w, String description){
		leftDownPanelFunctions.add(w,description,true);
	}
}