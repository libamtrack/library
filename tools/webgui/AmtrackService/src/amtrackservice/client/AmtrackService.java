package amtrackservice.client;

import java.util.ArrayList;
import java.util.HashMap;

import amtrackservice.client.gui.MainView;
import amtrackservice.client.parser.ConfigParser;
import amtrackservice.client.parser.TemplateParser;
import amtrackservice.shared.Config;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.visualization.client.VisualizationUtils;
import com.google.gwt.visualization.client.visualizations.PieChart;
import com.google.gwt.xml.client.impl.DOMParseException;

/**
 * This is the entry point for the webservice. Entry point classes define
 * <code>onModuleLoad()</code>.
 * 
 * @author Christoph Kolb, 2010 Hochschule Heilbronn
 */
public class AmtrackService implements EntryPoint {

	private MainView gui;
	private Config config = null;
	private static final String URL_CONFIG = "config.xml";
	private static final String URL_TEMPLATES = "templates/";

	private ArrayList<Calculation> openCalculations = new ArrayList<Calculation>();

	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
		gui = new MainView(this);
		initConfig();
		Logger.info("Service ready...");

//		Runnable onLoadCallback = new Runnable() {
//			public void run() {}
//		};
//		VisualizationUtils.loadVisualizationApi(onLoadCallback, PieChart.PACKAGE);
	}

	/**
	 * opens the calculation form defined by the name
	 * 
	 * @param method
	 *            name of the calculation form
	 */
	public void openCalculation(String method) {
		Logger.info("opening " + method);
		fetchTemplate(method);
	}

	/**
	 * this method is invoked by ServiceControl, when the template has been
	 * parsed to a Calculation Object. a tab will be added to the gui.
	 * 
	 * @param calc
	 */
	public void openCalculation(Calculation calc) {
		openCalculations.add(calc);
		gui.addTabPanel(calc.getCalculationPanel(), calc.getName());
		Logger.info(calc.getName() + " has been loaded");
	}

	/**
	 * initializes the gui. the list of functions in the menu will be filled.
	 */
	public void initGui() {
		gui.initMenu(getConfig().getTemplates().keySet().toArray(
				new String[0]));
	}

	
	public void fetchTemplate(String templateName) {
		String file = GWT.getModuleBaseURL() + URL_TEMPLATES + config.getTemplates().get(templateName);
		Logger.info("Requesting file: "+file+" from Server");
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET,
				file);

		try {
			requestBuilder.sendRequest(null, new RequestCallback() {
				public void onError(Request request, Throwable exception) {
					requestFailed(exception);
				}

				public void onResponseReceived(Request request,
						Response response) {
					try {
						Calculation calc = TemplateParser.readTemplate(response.getText());
						openCalculation(calc);
					} catch (DOMParseException e) {
						Logger.error("Couldn't parse template file");
						Logger.log(e.getMessage());
					}
				}
			});
		} catch (RequestException ex) {
			requestFailed(ex);
		} catch (DOMParseException e) {
			parseFailed(e);
		}
	}

	public void initConfig() {		
		String file = GWT.getModuleBaseURL() + URL_CONFIG;
		Logger.info("Requesting file: "+file+" from Server");
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET,
				file);
		try {
			requestBuilder.sendRequest(null, new RequestCallback() {
				public void onError(Request request, Throwable exception) {
					requestFailed(exception);
				}

				public void onResponseReceived(Request request,
						Response response) {
					try {
						Config config = ConfigParser.readConfig(response.getText());
						setConfig(config);
						if(config != null){
							initGui();
						}
						CalculationControl.getInstance().setConfig(config);
					} catch (DOMParseException e) {
						Logger.error("Couldn't parse config file");
						Logger.log(e.getMessage());
					}
				}
			});
		} catch (RequestException ex) {
			requestFailed(ex);
		} catch (DOMParseException e) {
			parseFailed(e);
		}
	}
	private void requestFailed(Throwable exception) {
		Logger.error("Failed to retrieve " + URL_CONFIG
				+ exception.getMessage());
	}

	private void parseFailed(Throwable exception) {
		Logger.error("Failed to parse: " + exception.getMessage());
	}

	public HashMap<String, String> getCalculationResult(long CalculationId) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public void setConfig(Config config){
		this.config = config;
	}

	public Config getConfig() {
		return config;
	}

	public void fetchCalculation(String templateName, long id) {
		String file = GWT.getModuleBaseURL() + URL_TEMPLATES + config.getTemplates().get(templateName);
		Logger.info("Requesting file: "+file+" from Server");
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET,
				file);

		try {
			requestBuilder.sendRequest(null, new RequestCallback() {
				public void onError(Request request, Throwable exception) {
					requestFailed(exception);
				}

				public void onResponseReceived(Request request,
						Response response) {
					try {
						Calculation calc = TemplateParser.readTemplate(response.getText());
						CalculationControl.getInstance().getCalculationResult(calc);
						openCalculation(calc);
					} catch (DOMParseException e) {
						Logger.error("Couldn't parse template file");
						Logger.log(e.getMessage());
					}
				}
			});
		} catch (RequestException ex) {
			requestFailed(ex);
		} catch (DOMParseException e) {
			parseFailed(e);
		}		
	}
}
