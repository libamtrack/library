package amtrackservice.client.parser;

import java.util.HashMap;

import amtrackservice.shared.Config;

import com.google.gwt.xml.client.Element;
import com.google.gwt.xml.client.Node;
import com.google.gwt.xml.client.NodeList;
import com.google.gwt.xml.client.impl.DOMParseException;

public class ConfigParser extends AbstractXMLParser {

	private static final String KEYWORD_PATH = "path";
	private static final String KEYWORD_TEMPLATE = "template";
	private static final String KEYWORD_NAME = "name";
	private static final String KEYWORD_CONFIG = "config";
	private static final String KEYWORD_FILE = "file";

	public static Config readConfig(String xmlTree) throws DOMParseException {
		Element root = parseXml(xmlTree);
		if (!root.getNodeName().equals(KEYWORD_CONFIG))
			throw new DOMParseException("Seems not to be a valid Config file");
		NodeList childNodes = root.getElementsByTagName(KEYWORD_TEMPLATE);
		HashMap<String, String> templates = new HashMap<String, String>();
		Node n = null;
		for (int i = 0; i < childNodes.getLength(); i++) {
			n = childNodes.item(i);
			HashMap<String, String> attMap = readAttributes(n);
			templates.put(attMap.get(KEYWORD_NAME), attMap.get(KEYWORD_FILE));
		}
		return new Config(readAttributes(root).get(KEYWORD_PATH), templates);

	}

}
