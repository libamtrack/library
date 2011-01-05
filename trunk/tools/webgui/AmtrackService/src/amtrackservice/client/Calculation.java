package amtrackservice.client;

import java.util.ArrayList;
import java.util.HashMap;

import amtrackservice.client.AmtrackService.AmtrackServiceResources;
import amtrackservice.client.gui.elements.AmWidget;

import com.google.gwt.core.client.GWT;
import com.google.gwt.dom.client.Style.Unit;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.event.dom.client.KeyCodes;
import com.google.gwt.event.dom.client.KeyPressEvent;
import com.google.gwt.event.dom.client.KeyPressHandler;
import com.google.gwt.event.logical.shared.CloseEvent;
import com.google.gwt.event.logical.shared.CloseHandler;
import com.google.gwt.user.client.Timer;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.DecoratedPopupPanel;
import com.google.gwt.user.client.ui.DialogBox;
import com.google.gwt.user.client.ui.DockLayoutPanel;
import com.google.gwt.user.client.ui.Grid;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.HasVerticalAlignment;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.PopupPanel;
import com.google.gwt.user.client.ui.PushButton;
import com.google.gwt.user.client.ui.ScrollPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.Widget;

/**
 * Represents a Calculation and generates the input/output gui
 * @author Christoph Kolb, 2010 Hochschule Heilbronn
 */
public class Calculation {

	private String name;
	private String wrapper;
	private long timeID;
	private boolean calculated;
	private HTML description;
	
	private DecoratedPopupPanel loggerPopUp;
	private String loggerText;
	private ScrollPanel dialogContents;

	private ArrayList<AmWidget> inputWidgets = new ArrayList<AmWidget>();
	private ArrayList<AmWidget> outputWidgets = new ArrayList<AmWidget>();

	/**
	 * 
	 * @return List of output widgets
	 */
	public ArrayList<AmWidget> getInputWidgets() {
		return inputWidgets;
	}

	/**
	 * 
	 * @return List of input widgets
	 */
	public ArrayList<AmWidget> getOutputWidgets() {
		return outputWidgets;
	}

	/**
	 * 
	 * @param timeID ID
	 * @param name will be displayed in the title
	 * @param wrapper the name of the function wrapper
	 * @param description a description of this calculation
	 * @param inputWidgets
	 * @param outputWidgets
	 * @param calculated
	 */
	public Calculation(long timeID, String name, String wrapper,
			HTML description, ArrayList<AmWidget> inputWidgets,
			ArrayList<AmWidget> outputWidgets, boolean calculated) {
		this.name = name;
		this.wrapper = wrapper;
		this.description = description;
		this.inputWidgets = inputWidgets;
		this.outputWidgets = outputWidgets;
		this.calculated = calculated;
		this.timeID = timeID;
		Logger.init(this);
		
		loggerPopUp = new DecoratedPopupPanel();
		loggerPopUp.setGlassEnabled(true);
		loggerPopUp.setAnimationEnabled(true);
		loggerPopUp.setAutoHideEnabled(true);
		loggerPopUp.setWidth("720px");
		
	}

	/**
	 * 
	 * @param name
	 * @param wrapper
	 * @param description
	 * @param inputWidgets
	 * @param outputWidgets
	 */
	public Calculation(String name, String wrapper, HTML description,
			ArrayList<AmWidget> inputWidgets, ArrayList<AmWidget> outputWidgets) {
		this(0, name, wrapper, description, inputWidgets, outputWidgets, false);
	}

	/**
	 * 
	 * @return the Panel, containing the input/output form
	 */
	public DockLayoutPanel getCalculationPanel() {
		Logger.init(this);
		
		AmtrackServiceResources resources = GWT.create(AmtrackServiceResources.class);
		
		DockLayoutPanel panel = new DockLayoutPanel(Unit.MM);
		HorizontalPanel forms = new HorizontalPanel();
		VerticalPanel head = new VerticalPanel();
		Grid input = new Grid(inputWidgets.size() + 4, 2);
		Grid output = new Grid(outputWidgets.size() + 1, 1);

		Image loadDefaultsImage = new Image(resources.defaults());		
		loadDefaultsImage.setSize("48px", "48px");
		PushButton defaultsButton = new PushButton(loadDefaultsImage);
		defaultsButton.setHeight("48px");
		defaultsButton.setWidth("48px");
		defaultsButton.addClickHandler(new ClickHandler() {
			public void onClick(ClickEvent event) {
				loadDefaults();
			}
		});

		Image reloadImage = new Image(resources.reload());
		reloadImage.setSize("48px", "48px");
		PushButton reloadButton = new PushButton(reloadImage);
		reloadButton.setHeight("48px");
		reloadButton.setWidth("48px");
		reloadButton.addClickHandler(new ClickHandler() {
			public void onClick(ClickEvent event) {
				recalculate();
			}
		});

		Image addImage = new Image(resources.add());
		addImage.setSize("48px", "48px");
		PushButton addButton = new PushButton(addImage);
		addButton.setHeight("48px");
		addButton.setWidth("48px");
		addButton.addClickHandler(new ClickHandler() {
			public void onClick(ClickEvent event) {
				//recalculate();
			}
		});

		Image statusImage = new Image(resources.status());
		statusImage.setSize("48px", "48px");
		PushButton statusButton = new PushButton(statusImage);
		statusButton.setHeight("48px");
		statusButton.setWidth("48px");
		statusButton.addClickHandler(new ClickHandler() {
			public void onClick(ClickEvent event) {
				showStatus();
			}
		});
		
		statusButton.addKeyPressHandler(new KeyPressHandler() {
			
			@Override
			public void onKeyPress(KeyPressEvent event) {
				if( event.getCharCode() == KeyCodes.KEY_ESCAPE ){
					hideStatus();
				}				
			}
		});
		
		head.add(description);

		ScrollPanel scroll = new ScrollPanel();
		forms.add(input);
		forms.add(output);
		scroll.add(forms);

		panel.addNorth(head, 5);
		panel.add(scroll);

		input.setCellPadding(0);
		input.setCellSpacing(0);

		output.setCellPadding(0);
		output.setCellSpacing(0);

		int inputRow = 0;
		for (AmWidget widget : inputWidgets) {
			if( widget.getClass().getName() == "amtrackservice.client.gui.elements.AmCombo"){		
				widget.getWidget().setWidth("150px");
			} else if (widget.getClass().getName() == "amtrackservice.client.gui.elements.AmTextField"){
				widget.getWidget().setWidth("150px");
			} 
			HorizontalPanel desc = new HorizontalPanel();
			desc.setVerticalAlignment(HasVerticalAlignment.ALIGN_MIDDLE);
			desc.add(new HTML("<p align=\"right\">" + widget.getLabel().getText() + "&nbsp;&nbsp;&nbsp;</p>"));
			
			Image helpImage = new Image(resources.help());
			helpImage.setSize("16px", "16px");
			PushButton helpButton = new PushButton(helpImage);
			helpButton.setHeight("16px");
			helpButton.setWidth("16px");			
			final String helpText = widget.getDescription().toString();
			helpButton.addClickHandler(new ClickHandler() {
				public void onClick(ClickEvent event) {
			          Window.confirm(helpText);
				}
			});
			//desc.add(helpButton);			
			input.setWidget(inputRow, 0, desc);
			input.setWidget(inputRow, 1, widget.getWidget());
			
			input.getCellFormatter().setHorizontalAlignment(inputRow, 0, HasHorizontalAlignment.ALIGN_RIGHT);
			
			inputRow++;
		}
		input.setWidget(inputRow, 0, new HTML("<p align=\"right\"><b>RECALCULATE</b>&nbsp;&nbsp;&nbsp;</p>"));
		input.setWidget(inputRow, 1, reloadButton);
		
//		input.setWidget(inputRow+1, 0, new HTML("<p align=\"right\"><b>ADD TO PLOT</b>&nbsp;&nbsp;&nbsp;</p>"));
//		input.setWidget(inputRow+1, 1, addButton);
		
		input.setWidget(inputRow+1, 0, new HTML("<p align=\"right\"><b>LOAD DEFAULTS</b>&nbsp;&nbsp;&nbsp;</p>"));
		input.setWidget(inputRow+1, 1, defaultsButton);
		
		input.setWidget(inputRow+2, 0, new HTML("<p align=\"right\"><b>STATUS</b>&nbsp;&nbsp;&nbsp;</p>"));
		input.setWidget(inputRow+2, 1, statusButton);
		
		input.getCellFormatter().setHeight(inputRow, 1, "65px");
		input.getCellFormatter().setHeight(inputRow+1, 1, "65px");
		input.getCellFormatter().setHeight(inputRow+2, 1, "65px");
		
		int outputRow = 1;
		for (AmWidget widget : outputWidgets) {
			output.setWidget(outputRow, 0, widget.getWidget());
			outputRow++;
		}

	    // Create a table to layout the content
	    this.dialogContents = new ScrollPanel();
	    this.dialogContents.setHeight("500px");
	    this.loggerPopUp.setWidget(this.dialogContents);
		
		return panel;
	}
	
	public void recalculate(){
	    this.clearLoggerMessages();
		calculate();
	    // repeat until Calculation is done or n attempts
		// were performed
	    Timer t = new Timer() {
	    	int n = 40;	    	
	    	public void run() {
	    		n--;
	    		if( isCalculated() || n < 0 ){
	    			this.cancel();  // cancel timer
	    		} else {
	    			refresh(); // refresh calculation
	    		}
	        }
	    };
	    // repeating interval is specified below
	    t.scheduleRepeating(250);
	}

	private void loadDefaults() {
		for (AmWidget widget : inputWidgets) {
			widget.setDefault();
		}
	}
	
	private void calculate() {
		CalculationControl.getInstance().calculate(this);
	}

	public void showStatus() {
		this.dialogContents.setWidget(new HTML(this.loggerText));
		this.loggerPopUp.center();
		this.loggerPopUp.show();
	}	
	
	public void hideStatus() {
		this.loggerPopUp.hide();
	}
	
	/**
	 * sets an id to this calculation
	 * @param timeID
	 */
	public void setTimeID(long timeID) {
		this.timeID = timeID;
	}

	private void refresh() {
		CalculationControl.getInstance().getCalculationResult(this);
	}

	/**
	 * fills all output widgets of this calculation with the result
	 * @param result
	 */
	public void setCalculationResult(HashMap<String, String> result) {
		if (result != null) {
			for (AmWidget widget : outputWidgets) {
				widget.setValue(result);
			}
		}
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getWrapper() {
		return wrapper;
	}

	public void setWrapper(String wrapper) {
		this.wrapper = wrapper;
	}

	public long getTimeID() {
		return timeID;
	}

	public boolean isCalculated() {
		return calculated;
	}

	public void setCalculated(boolean calculated) {
		this.calculated = calculated;
	}

	public HTML getDescription() {
		return description;
	}

	public void setDescription(HTML description) {
		this.description = description;
	}
	
	public void clearLoggerMessages(){
		this.loggerText = "<H1><B>Click outside this window to displace it</B></H1>";
	}
	
	public void addLoggerMessage(String message){
		this.loggerText += "<BR>" + message;
	}
}
