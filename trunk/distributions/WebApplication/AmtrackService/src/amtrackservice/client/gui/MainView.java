package amtrackservice.client.gui;

import amtrackservice.client.AmtrackService;
import amtrackservice.client.AmtrackService.AmtrackServiceResources;
import amtrackservice.client.BuildInfo;

import com.google.gwt.core.client.GWT;
import com.google.gwt.dom.client.Style.Unit;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.DOM;
import com.google.gwt.user.client.ui.DecoratedStackPanel;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.HasVerticalAlignment;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.Label;
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
	private HTML rightMainTab;
	
	/**
	 * Constructor of MainView
	 * @param service 
	 */
	public MainView(AmtrackService service) {
		
		AmtrackServiceResources resources = GWT.create(AmtrackServiceResources.class);
		
		/**  Constructors **/
		RootLayoutPanel rootPanel  = RootLayoutPanel.get();
	    HorizontalPanel wholePagePanel = new HorizontalPanel();
		VerticalPanel leftPanel = new VerticalPanel();
	    VerticalPanel rightPanel = new VerticalPanel();

	    Image leftPanelUpperLogo = new Image(resources.logo());
		leftDownPanelFunctions = new DecoratedStackPanel();
	    
		rightMainTab = new HTML("<center><big><font size=\"+2\"><br><br>Welcome to <a href=\"http://libamtrack.dkfz.org\">libamtrack</a> webGUI</font></big><br><br>" + 
					"<big>HOWTO</big><br><br></center>" +
					"<p style=\"margin-left:6cm\"><font size=\"2\">1. Choose a <B>function</B> from the left panel.<br>"+
					"2. Enter desired values for the <B>input parameters</B> or use defaults.<br>"+
					"3. Press <B>\"Recalculate\"</B> to start new computation (be careful: this will erase all previous computations).<br>"+
					"4. Press <B>\"Add plot\"</B> to add a new computation to your plot (note that every computation needs a unique label).<br>"+
					"5. Press <B>\"Load default\"</B> to reset input values.<br>"+
					"6. Press <B>\"Status\"</B> to get information on computation on the server.<br>"+
					"7. Results will be printed <B>on the right.</B><br>"+
					"8. Use the <B>questionmark buttons</B> to get more information on the input parameters.<br>"+
					"9. If you experience problems, please try to clean your <B>browser cache</B> first. If problems persist, please contact:<br>"+
		            "<b><u><a href=\"mailto:Leszek.Grzanka@ifj.edu.pl\" target=\"_blank\">Leszek.Grzanka@ifj.edu.pl</a></u> or <u><a href=\"mailto:Steffen.Greilich@dkfz.de\" target=\"_blank\">Steffen.Greilich@dkfz.de</a></u></b><br></font>" +
		            "</p><center><I>libamtrack version: " + BuildInfo.revisionNumber + "</I>" + 
		            "</center>");
		
	    centralTabPanel = new TabLayoutPanel(32, Unit.PX);

	    
		/**  Tree of life **/
		rootPanel.add(wholePagePanel);		
		wholePagePanel.add(leftPanel);		
	    wholePagePanel.add(rightPanel);	    
		
	    leftPanel.add(leftPanelUpperLogo);
	    leftPanel.add(leftDownPanelFunctions);

	    rightPanel.add(centralTabPanel);
	    
		
		/**  Dimensions **/
	    rootPanel.setSize("1180px", "800px");
		rootPanel.setWidgetLeftWidth(wholePagePanel, 0.0, Unit.PX, 1180.0, Unit.PX);
		rootPanel.setWidgetTopHeight(wholePagePanel, 0.0, Unit.PX, 800.0, Unit.PX);

	    rightPanel.setWidth("600px");
		
		leftPanel.setWidth("180px");
	    leftPanel.setCellHeight(leftPanelUpperLogo, "120px");
	    leftPanel.setCellWidth(leftDownPanelFunctions, "180px");
		
	    leftPanelUpperLogo.setSize("180px", "120px");
	    leftDownPanelFunctions.setWidth("180px");	    

	    centralTabPanel.setSize("1000px", "800px");
		
		
		/**  Alignment **/
		wholePagePanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_CENTER);

		leftPanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_CENTER);
	    leftPanel.setCellVerticalAlignment(leftPanelUpperLogo, HasVerticalAlignment.ALIGN_BOTTOM);
	    leftPanel.setCellHorizontalAlignment(leftPanelUpperLogo, HasHorizontalAlignment.ALIGN_CENTER);
	    leftPanel.setCellHorizontalAlignment(leftDownPanelFunctions, HasHorizontalAlignment.ALIGN_CENTER);
	    
	    rightPanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_LEFT);
		
	    
		/**  Others **/
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
		
		centralTabPanel.add(input, newLabel(input, name));
		centralTabPanel.selectTab(input);		
	}
	
	/**
	 * TODO
	 * @param widget
	 * @param string
	 * @return
	 */
	private Widget newLabel(final Widget widget, final String string) {
		final HorizontalPanel tabLabelHPanel = new HorizontalPanel();
		final Label label = new Label(string);
		DOM.setStyleAttribute(label.getElement(), "whiteSpace", "nowrap");

		AmtrackServiceResources resources = GWT.create(AmtrackServiceResources.class);
		PushButton closeButton = new PushButton(new Image(resources.close()));
		closeButton.setHeight("16px");
		closeButton.setWidth("16px");
		closeButton.addClickHandler(new ClickHandler() {
			public void onClick(ClickEvent event) {
				int index = centralTabPanel.getSelectedIndex();
				if( index > -1 ){
					centralTabPanel.remove(index);
				}
				if( centralTabPanel.getWidgetCount() == 0){
					centralTabPanel.add(rightMainTab);
				}
			}
		});
		tabLabelHPanel.add(label);
		tabLabelHPanel.setCellVerticalAlignment(label, HasVerticalAlignment.ALIGN_MIDDLE);
		tabLabelHPanel.add(new HTML("&nbsp&nbsp"));
		tabLabelHPanel.add(closeButton);
		return tabLabelHPanel;
	} 
	
	/**
	 * TODO
	 * @param widget
	 * @param description
	 */
	public void addWidgetToLeftPanel(Widget widget, String description){
		leftDownPanelFunctions.add(widget,description,true);
	}

}