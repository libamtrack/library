package amtrackservice.client;

import java.util.ArrayList;
import java.util.HashMap;

import amtrackservice.client.gui.MainView;
import amtrackservice.client.parser.ConfigParser;
import amtrackservice.client.parser.TemplateParser;
import amtrackservice.shared.Config;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.GWT;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.http.client.Request;
import com.google.gwt.http.client.RequestBuilder;
import com.google.gwt.http.client.RequestCallback;
import com.google.gwt.http.client.RequestException;
import com.google.gwt.http.client.Response;
import com.google.gwt.resources.client.ClientBundle;
import com.google.gwt.resources.client.ImageResource;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.HasVerticalAlignment;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.Label;
import com.google.gwt.user.client.ui.PushButton;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.visualization.client.VisualizationUtils;
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

	public interface AmtrackServiceResources extends ClientBundle {
		  @Source("Play.png")
		  ImageResource play();
		  
		  @Source("Close.png")
		  ImageResource close();		  
		  
		  @Source("Logo.png")
		  ImageResource logo();
		  
		  @Source("Reload.png")
		  ImageResource reload();
		  
		  @Source("Defaults.png")
		  ImageResource defaults();

		  @Source("Add.png")
		  ImageResource add();

		  @Source("Help.png")
		  ImageResource help();
		  
		  @Source("Status.png")
		  ImageResource status();
	}
	
	/**
	 * This is the entry point method.
	 */
	public void onModuleLoad() {
		gui = new MainView(this);
		initConfig();
		Logger.info("Service ready...");

		Runnable onLoadCallback = new Runnable() {
			public void run() {}
		};
		
		VisualizationUtils.loadVisualizationApi(onLoadCallback, "corechart");
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
		calc.clearLoggerMessages();
		openCalculations.add(calc);
		gui.addTabPanel(calc.getCalculationPanel(), calc.getName());
		Logger.info(calc.getName() + " has been loaded");
		calc.recalculate();
	}

	/**
	 * initializes the gui. the list of functions in the menu will be filled.
	 */
	public void initGui() {
		/** TODO move this method to GUI class */
		AmtrackServiceResources resources = GWT.create(AmtrackServiceResources.class);
		
		VerticalPanel generalFunctionPanel = new VerticalPanel();
		generalFunctionPanel.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_LEFT);
		generalFunctionPanel.setSpacing(5);
		for(String m : getConfig().getTemplates().keySet().toArray(new String[0])){
				HorizontalPanel openCalculationPanel = new HorizontalPanel();
			    openCalculationPanel.setSpacing(3);
			    openCalculationPanel.setVerticalAlignment(HasVerticalAlignment.ALIGN_MIDDLE);
			    final String ms = m;
			    ClickHandler ch = new ClickHandler() {					
					@Override
					public void onClick(ClickEvent event) {
						openCalculation(ms);
					}
				};
			    PushButton openCalculationButton = new PushButton(new Image(resources.play()),ch);
			    openCalculationPanel.add(openCalculationButton);
			    HTML openCalculationText = new HTML(m);
			    openCalculationPanel.add(openCalculationText);
			    generalFunctionPanel.add(openCalculationPanel);
		}			
		gui.addWidgetToLeftPanel(generalFunctionPanel,"General calculations");
		gui.addWidgetToLeftPanel(new Label("none yet"),"Others");
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
		RequestBuilder requestBuilder = new RequestBuilder(RequestBuilder.GET, file);

		try {
			requestBuilder.sendRequest(null, new RequestCallback() {
				public void onError(Request request, Throwable exception) {
					requestFailed(exception);
				}

				public void onResponseReceived(Request request, Response response) {
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
