package amtrackservice.client.parser;

import java.util.ArrayList;
import java.util.HashMap;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.xml.client.Element;
import com.google.gwt.xml.client.Node;
import com.google.gwt.xml.client.NodeList;
import com.google.gwt.xml.client.impl.DOMParseException;

import amtrackservice.client.Calculation;
import amtrackservice.client.MapList;
import amtrackservice.client.gui.elements.AmCombo;
import amtrackservice.client.gui.elements.AmInputList;
import amtrackservice.client.gui.elements.AmList;
import amtrackservice.client.gui.elements.AmPlot;
import amtrackservice.client.gui.elements.AmTextField;
import amtrackservice.client.gui.elements.AmWidget;

public class TemplateParser extends AbstractXMLParser {

	private static final String KEYWORD_DESCRIPTION = "description";
	private static final String KEYWORD_INPUT = "input";
	private static final String KEYWORD_OUTPUT = "output";
	private static final String KEYWORD_ROOT = "root";
	private static final String KEYWORD_NAME = "name";
	private static final String KEYWORD_WRAPPER = "wrapper";
	private static final String KEYWORD_LABEL = "label";
	private static final String KEYWORD_DATATYPE = "datatype";
	private static final String KEYWORD_TYPE = "type";
	private static final String KEYWORD_PRESET = "preset";
	private static final String KEYWORD_ITEM = "item";
	private static final String KEYWORD_DATA = "data";
	private static final String KEYWORD_ENTRY = "entry";
	private static final String KEYWORD_KEY = "key";
	private static final String KEYWORD_VALUE = "value";
	private static final String KEYWORD_ELEMENT = "element";
	private static final String KEYWORD_AXIS = "axis";
	private static final String KEYWORD_X = "x";
	private static final String KEYWORD_Y = "y";
	private static final String KEYWORD_Z = "z";

	private static final String GUI_INPUTLIST = "inputlist";
	private static final String GUI_LIST = "list";
	private static final String GUI_COMBO = "combo";
	private static final String GUI_TEXTFIELD = "field";
	private static final String GUI_PLOT = "plot";

	public static Calculation readTemplate(String xmlTree)
			throws DOMParseException {
		Element root = parseXml(xmlTree);
		if (!root.getNodeName().equals(KEYWORD_ROOT))
			throw new DOMParseException("Seems not to be a valid Template file");

		String name = "";
		String wrapper = "";
		HTML description = new HTML();
		ArrayList<AmWidget> inputWidgets = new ArrayList<AmWidget>();
		ArrayList<AmWidget> outputWidgets = new ArrayList<AmWidget>();

		HashMap<String, String> attributes = readAttributes(root);
		name = attributes.get(KEYWORD_NAME);
		wrapper = attributes.get(KEYWORD_WRAPPER);

		NodeList childNodes = root.getChildNodes();
		Node n = null;
		for (int i = 0; i < childNodes.getLength(); i++) {
			n = childNodes.item(i);
			if (n.getNodeName().equals(KEYWORD_DESCRIPTION)) {
				description = readDescription(n);
			} else if (n.getNodeName().equals(KEYWORD_INPUT)) {
				inputWidgets = readGuiElements(n);
			} else if (n.getNodeName().equals(KEYWORD_OUTPUT)) {
				outputWidgets = readGuiElements(n);
			}
		}
		return new Calculation(name, wrapper, description, inputWidgets,
				outputWidgets);
	}

	private static HTML readDescription(Node descriptionNode) {
		return translateHTML(getTextValue(descriptionNode));
	}

	private static ArrayList<AmWidget> readGuiElements(Node elementNodes) {
		ArrayList<AmWidget> widList = new ArrayList<AmWidget>();
		NodeList childNodes = elementNodes.getChildNodes();
		Node n = null;
		for (int i = 0; i < childNodes.getLength(); i++) {
			n=childNodes.item(i);
			if (n.getNodeName().equals(KEYWORD_ELEMENT)) {
				AmWidget element = createGuiElement(n);
				if (element != null)
					widList.add(element);
			}
		}
		return widList;
	}

	private static HTML translateHTML(String htmlText) {
		return new HTML(htmlText.replace('[', '<').replace(']', '>'));
	}

	private static AmWidget createGuiElement(Node guiElement) {
		AmWidget element = null;
		HashMap<String, String> attributes = readAttributes(guiElement);
		String label = attributes.get(KEYWORD_LABEL);
		String type = attributes.get(KEYWORD_TYPE);
		String datatype = attributes.get(KEYWORD_DATATYPE);

		HTML description = null;
		MapList<String, String> entry = null;
		MapList<String, String> preset = null;
		String dataX = null;
		String dataY = null;
		String dataZ = null;
		
		boolean xAxisLog = false;
		boolean yAxisLog = false;
		
		NodeList childNodes = guiElement.getChildNodes();
		Node n = null;
		for (int i = 0; i < childNodes.getLength(); i++) {
			n=childNodes.item(i);
			if (n.getNodeName().equals(KEYWORD_DESCRIPTION)) {
				description = readDescription(n);
			} else if (n.getNodeName().equals(KEYWORD_PRESET)) {
				preset = readItemList(n);
			} else if (n.getNodeName().equals(KEYWORD_DATA)) {
				HashMap<String, String> data_attributes = readAttributes(n);
				dataX = data_attributes.get(KEYWORD_X);
				dataY = data_attributes.get(KEYWORD_Y);
				dataZ = data_attributes.get(KEYWORD_Z);
			} else if (n.getNodeName().equals(KEYWORD_ENTRY)) {
				entry = readItemList(n);
			} else if (n.getNodeName().equals(KEYWORD_AXIS)) {
				HashMap<String, String> data_attributes = readAttributes(n);
				String xAxisType = data_attributes.get(KEYWORD_X);
				if( (xAxisType != null) && (xAxisType.equalsIgnoreCase("log")) ){
					xAxisLog = true;
				}
				String yAxisType = data_attributes.get(KEYWORD_Y);
				if( (yAxisType != null) && (yAxisType.equalsIgnoreCase("log")) ){
					yAxisLog = true;
				}
			}
		}

		if (type.equals(GUI_COMBO)) {
			element = new AmCombo(label, datatype, description, entry, preset,
					dataX);
		} else if (type.equals(GUI_TEXTFIELD)) {
			element = new AmTextField(label, datatype, description, preset,
					dataX);
		} else if (type.equals(GUI_INPUTLIST)) {
			element = new AmInputList(label, datatype, description, preset,
					dataX);
		} else if (type.equals(GUI_LIST)) {
			element = new AmList(label, datatype, description, preset, dataX);
		} else if (type.equals(GUI_PLOT)) {
			element = new AmPlot(label, datatype, description, preset, dataX, dataY, dataZ, xAxisLog, yAxisLog);
		}
		return element;
	}

	private static MapList<String, String> readItemList(Node itemListNode) {
		MapList<String, String> itemList = new MapList<String, String>();
		NodeList childNodes = itemListNode.getChildNodes();
		Node n = null;
		for (int i = 0; i < childNodes.getLength(); i++) {
			n = childNodes.item(i);
			HashMap<String, String> attributes = readAttributes(n);
			if (attributes.get(KEYWORD_KEY) != null)
				itemList.put(attributes.get(KEYWORD_KEY), attributes
						.get(KEYWORD_VALUE));
			else
				itemList.put(attributes.get(KEYWORD_VALUE), attributes
						.get(KEYWORD_VALUE));
		}
		return itemList;
	}
}
