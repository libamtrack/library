package amtrackservice.shared;

import java.io.Serializable;
import java.util.HashMap;

import amtrackservice.client.parser.AbstractXMLParser;

public class Config extends AbstractXMLParser implements Serializable{

	private static final long serialVersionUID = 928505004033338560L;
	private HashMap<String, String> templateFilenames = new HashMap<String, String>();
	private HashMap<String, String> groupNames = new HashMap<String, String>();
	private String serverPath;
	private static final String PATH_WRAPPER = "wrapper/";
	private static final String PATH_RESULTS = "result/";

	public Config(String serverPath, HashMap<String, String> templateFilenames, HashMap<String, String> groupNames) {
		this.serverPath = serverPath;
		this.templateFilenames = templateFilenames;
		this.groupNames = groupNames;
	}

	public HashMap<String, String> getTemplateFilenames() {
		return templateFilenames;
	}

	public HashMap<String, String> getGroupNames() {
		return groupNames;
	}

	public String getWrapperPath() {
		return serverPath + PATH_WRAPPER;
	}
	
	public String getResultPath() {
		return serverPath + PATH_RESULTS;
	}

}
