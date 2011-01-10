package amtrackservice.client.parser;

import java.util.HashMap;

import com.google.gwt.xml.client.Element;
import com.google.gwt.xml.client.NamedNodeMap;
import com.google.gwt.xml.client.Node;
import com.google.gwt.xml.client.XMLParser;
import com.google.gwt.xml.client.impl.DOMParseException;

public class AbstractXMLParser {

	protected static Element parseXml(String xmlText) throws DOMParseException{
		Element root = XMLParser.parse(xmlText).getDocumentElement();
		XMLParser.removeWhitespace(root);
		return root;
	}
	
	protected static String getTextValue(Node node){
		return node.getFirstChild().getNodeValue().trim();
	}
	
	protected static HashMap<String, String> readAttributes(Node xmlNode) {
		HashMap<String, String> attMap = new HashMap<String, String>();
		NamedNodeMap atts = xmlNode.getAttributes();
		for (int i = 0; i < atts.getLength(); i++) {
			Node att = atts.item(i);
			String name = att.getNodeName();
			attMap.put(name, att.getNodeValue().trim());
		}
		return attMap;
	}
}
