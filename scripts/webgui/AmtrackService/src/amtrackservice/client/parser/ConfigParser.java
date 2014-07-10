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
	private static final String KEYWORD_GROUP = "group";

	public static Config readConfig(String xmlTree) throws DOMParseException {
		Element root = parseXml(xmlTree);
		if (!root.getNodeName().equals(KEYWORD_CONFIG))
			throw new DOMParseException("Seems not to be a valid Config file");
		NodeList childNodes = root.getElementsByTagName(KEYWORD_TEMPLATE);
		HashMap<String, String> templateFilenames = new HashMap<String, String>();
		HashMap<String, String> groupNames = new HashMap<String, String>();
		Node n = null;
		for (int i = 0; i < childNodes.getLength(); i++) {
			n = childNodes.item(i);
			HashMap<String, String> attMap = readAttributes(n);
			templateFilenames.put(attMap.get(KEYWORD_NAME), attMap.get(KEYWORD_FILE));
			groupNames.put(attMap.get(KEYWORD_NAME), attMap.get(KEYWORD_GROUP));
		}
		return new Config(readAttributes(root).get(KEYWORD_PATH), templateFilenames, groupNames);

	}

}
