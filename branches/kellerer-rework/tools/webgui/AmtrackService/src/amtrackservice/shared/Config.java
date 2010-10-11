package amtrackservice.shared;

import java.io.Serializable;
import java.util.HashMap;

import amtrackservice.client.parser.AbstractXMLParser;

public class Config extends AbstractXMLParser implements Serializable{

	private static final long serialVersionUID = 928505004033338560L;
	private HashMap<String, String> templates = new HashMap<String, String>();
	private String serverPath;
	private static final String PATH_WRAPPER = "wrapper/";
	private static final String PATH_RESULTS = "result/";

	public Config(String serverPath, HashMap<String, String> templates) {
		this.serverPath = serverPath;
		this.templates = templates;
	}

	public HashMap<String, String> getTemplates() {
		return templates;
	}

	public String getWrapperPath() {
		return serverPath + PATH_WRAPPER;
	}
	
	public String getResultPath() {
		return serverPath + PATH_RESULTS;
	}

}
