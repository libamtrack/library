package amtrackservice.client;

import java.util.ArrayList;
import java.util.HashMap;

import amtrackservice.client.gui.elements.AmWidget;

import com.google.gwt.dom.client.Style.Unit;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.DockLayoutPanel;
import com.google.gwt.user.client.ui.Grid;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.ScrollPanel;
import com.google.gwt.user.client.ui.StackLayoutPanel;
import com.google.gwt.user.client.ui.VerticalPanel;

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

		DockLayoutPanel panel = new DockLayoutPanel(Unit.MM);
		HorizontalPanel forms = new HorizontalPanel();
		HorizontalPanel buttons = new HorizontalPanel();
		VerticalPanel head = new VerticalPanel();
		Grid input = new Grid(inputWidgets.size() + 1, 2);
		Grid output = new Grid(outputWidgets.size() + 1, 2);
		StackLayoutPanel descriptions = new StackLayoutPanel(Unit.MM);
		Button defaults = new Button("load defaults");
		defaults.addClickHandler(new ClickHandler() {
			@Override
			public void onClick(ClickEvent event) {
				loadDefaults();
			}
		});
		Button refresh = new Button("fetch result");
		refresh.addClickHandler(new ClickHandler() {
			@Override
			public void onClick(ClickEvent event) {
				refresh();

			}
		});
		Button calculate = new Button("start calculation");
		calculate.addClickHandler(new ClickHandler() {
			@Override
			public void onClick(ClickEvent event) {
				calculate();
			}
		});

		buttons.add(defaults);
		buttons.add(calculate);
		buttons.add(refresh);

		head.add(description);
		head.add(buttons);

		ScrollPanel scroll = new ScrollPanel();
		forms.add(input);
		forms.add(output);
		scroll.add(forms);

		panel.addWest(descriptions, 40);
		panel.addNorth(head, 30);
		panel.add(scroll);

		input.setWidget(0, 0, new HTML("<b>Input</b>"));
		output.setWidget(0, 0, new HTML("<b>Output</b>"));

		int inputRow = 1;
		for (AmWidget widget : inputWidgets) {
			input.setWidget(inputRow, 0, new HTML("<p align=\"right\">"
					+ widget.getLabel().getText() + "</p>"));
			input.setWidget(inputRow, 1, widget.getWidget());
			descriptions.add(widget.getDescription(), widget.getLabel()
					.getText(), 10);
			inputRow++;
		}

		int outputRow = 1;
		for (AmWidget widget : outputWidgets) {
			output.setWidget(outputRow, 0, new HTML("<p align=\"right\">"
					+ widget.getLabel().getText() + "</p>"));
			output.setWidget(outputRow, 1, widget.getWidget());
			descriptions.add(widget.getDescription(), widget.getLabel()
					.getText(), 10);
			outputRow++;
		}

		return panel;
	}

	private void loadDefaults() {
		for (AmWidget widget : inputWidgets) {
			widget.setDefault();
		}
	}
	
	private void calculate() {
		CalculationControl.getInstance().calculate(this);
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

}
