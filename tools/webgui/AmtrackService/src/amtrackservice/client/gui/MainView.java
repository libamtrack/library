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

	private TabLayoutPanel tabPanel;
	private DecoratedStackPanel leftPanel;

	/**
	 * Constructor of MainView
	 * @param service 
	 */
	public MainView(AmtrackService service) {
		Logger.init(this);
		
		AmtrackServiceResources resources = GWT.create(AmtrackServiceResources.class);
	    ClickHandler closeHandler = new ClickHandler() {					
			@Override
			public void onClick(ClickEvent event) {
				int index = tabPanel.getSelectedIndex();
				tabPanel.remove(index);
			}
		};
		
		HorizontalPanel horizontalPanel = new HorizontalPanel();
		horizontalPanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_CENTER);

		VerticalPanel westVerticalPanel = new VerticalPanel();
		westVerticalPanel.setWidth("300px");
		westVerticalPanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_CENTER);

		Image logoImage = new Image(resources.logo());
		logoImage.setSize("300px", "200px");
		
	    leftPanel = new DecoratedStackPanel();
	    leftPanel.setWidth("300px");

	    VerticalPanel verticalPanel_1 = new VerticalPanel();
	    verticalPanel_1.setSize("1000px", "");

	    HorizontalPanel horizontalPanel_1 = new HorizontalPanel();
	    horizontalPanel_1.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_RIGHT);
	    horizontalPanel_1.setSize("1000px", "32px");
	    
	    Image closeImage = new Image(resources.close());
	    closeImage.setSize("32px", "32px");

	    HTML htmlCloseCurrentCalculation = new HTML("Close current calculation: ", true);
	    htmlCloseCurrentCalculation.setStyleName("h1");
	    htmlCloseCurrentCalculation.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_RIGHT);
	    htmlCloseCurrentCalculation.setSize("200px", "32px");
	    
	    PushButton closeButton = new PushButton(closeImage,closeHandler);
	    closeButton.setSize("32px", "32px");

	    tabPanel = new TabLayoutPanel(10, Unit.MM);
	    tabPanel.setSize("1000px", "800px");
	    
		RootLayoutPanel rp  = RootLayoutPanel.get();
		rp.setSize("1300px", "1300px");
		
		rp.add(horizontalPanel);
		rp.setWidgetLeftWidth(horizontalPanel, 0.0, Unit.PX, 1800.0, Unit.PX);
		rp.setWidgetTopHeight(horizontalPanel, 0.0, Unit.PX, 1800.0, Unit.PX);
		
		horizontalPanel.add(westVerticalPanel);		
    
	    westVerticalPanel.add(logoImage);
	    westVerticalPanel.add(leftPanel);
	    westVerticalPanel.setCellHeight(logoImage, "200px");
	    westVerticalPanel.setCellVerticalAlignment(logoImage, HasVerticalAlignment.ALIGN_BOTTOM);
	    westVerticalPanel.setCellHorizontalAlignment(logoImage, HasHorizontalAlignment.ALIGN_CENTER);
	    westVerticalPanel.setCellWidth(leftPanel, "300px");
	    westVerticalPanel.setCellHorizontalAlignment(leftPanel, HasHorizontalAlignment.ALIGN_CENTER);

	    horizontalPanel.add(verticalPanel_1);	    
	    verticalPanel_1.add(horizontalPanel_1);	    
	    verticalPanel_1.setCellHorizontalAlignment(horizontalPanel_1, HasHorizontalAlignment.ALIGN_RIGHT);
	    	    
	    horizontalPanel_1.add(htmlCloseCurrentCalculation);
	    horizontalPanel_1.setCellHorizontalAlignment(htmlCloseCurrentCalculation, HasHorizontalAlignment.ALIGN_RIGHT);
	    horizontalPanel_1.setCellHeight(htmlCloseCurrentCalculation, "32px");
	    horizontalPanel_1.setCellVerticalAlignment(htmlCloseCurrentCalculation, HasVerticalAlignment.ALIGN_MIDDLE);
	    horizontalPanel_1.add(closeButton);
	    horizontalPanel_1.setCellHeight(closeButton, "32px");
	    horizontalPanel_1.setCellWidth(closeButton, "32px");
	    horizontalPanel_1.setCellHorizontalAlignment(closeButton, HasHorizontalAlignment.ALIGN_RIGHT);

	    verticalPanel_1.add(tabPanel);
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
	
	public void addWidgetToLeftPanel(Widget w, String description){
		leftPanel.add(w,description,true);
	}
}