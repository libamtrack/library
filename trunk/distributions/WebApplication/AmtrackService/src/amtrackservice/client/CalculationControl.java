package amtrackservice.client;

import java.util.HashMap;

import amtrackservice.client.gui.elements.AmWidget;
import amtrackservice.shared.Config;

import com.google.gwt.core.client.GWT;
import com.google.gwt.user.client.rpc.AsyncCallback;

/**
 * implements methods to contact the server.
 * @author Christoph Kolb, 2010 Hochschule Heilbronn
 *
 */
public class CalculationControl {

	/**
	 * Create a remote service proxy to talk to the server-side service.
	 */
	private final CalculationServiceAsync calculationService = GWT.create(CalculationService.class);

	/**
	 * Calculation configuration object
	 */
	private Config config = null;
	
	/**
	 * static instance
	 */
	private static CalculationControl instance = null;
	
	/**
	 * default constructor
	 */
	private CalculationControl() {	}

	/**
	 * getter for static instance
	 * @return static instance
	 */
	public static CalculationControl getInstance() {
		if (instance == null)
			instance = new CalculationControl();
		return instance;
	}

	/**
	 * TODO
	 * @param config
	 */
	public void setConfig(Config config) {
		this.config = config;
	}

	/**
	 * Request calculation results from server
	 * @param calculation
	 */
	public void getCalculationResult(final Calculation calculation) {
		Logger.init(calculation);
		Logger.info("Requesting result for calculation " + calculation.getName() + "(" + calculation.getTimeID() + ")");
		if( calculation.getTimeID() == 0 ){
			Logger.info("Process id not set yet...");
			//calculation.setCalculationResult(new HashMap<String, String>());
			calculation.setCalculated(false);
		} else {
			calculationService.requestResult(calculation.getTimeID(),
					new AsyncCallback<HashMap<String, String>>() {

				// TODO needs to be rewritten in clearer way
				public void onSuccess(HashMap<String, String> result) {						
					Logger.info("Reqest successful");
					if( result.get("exitcode") != null){
						String exitcode = result.get("exitcode").trim();
						if( exitcode == ""){
							Logger.info("Process not finished yet");
							//calculation.setCalculationResult(new HashMap<String, String>());
							calculation.setCalculated(false);
						}else{
							calculation.setCalculated(true);
							if( exitcode.equals("0")){
								Logger.info("Process finished, success" );
								String content = result.get("content");
								if( content == null)
									Logger.info("Output file content : null<br>");
								else{
									Logger.info("Output file content : <BR>" + content.replaceAll("(\r\n|\r|\n|\n\r)", "<br>"));
									int status = calculation.setCalculationResult(result);
									if( status != 0){
										Logger.error("Problem with filling output elements with data, code " + status);
										if( status == -1)
											Logger.error("Data X name in collection is different from data X name in serie");
										if( status == -2)
											Logger.error("Data Y name in collection is different from data Y name in serie");
										if( status == -3)
											Logger.error("Another data serie with the same name exists in collection");										
										if( status == -4)
											Logger.error("Couldn't prepare DataTable for plot");										
										Logger.show();
									}
								}
							} else {
								Logger.info("Standard output : <BR>" + result.get("stdout").replaceAll("(\r\n|\r|\n|\n\r)", "<br>"));
								Logger.info("Standard error : <BR>" + result.get("stderr").replaceAll("(\r\n|\r|\n|\n\r)", "<br>"));
								Logger.error("Process finished, failure, exit code = _" + exitcode + "_");
								int status = calculation.setCalculationResult(new HashMap<String, String>());
								if( status != 0){
									Logger.error("Problem with filling output elements with data");
									Logger.show();
								}
								calculation.showStatus();
							}
						}
					} else {
						calculation.setCalculated(false);
						Logger.info("Process not finished yet");
						//calculation.setCalculationResult(new HashMap<String, String>());
					}
				}

				public void onFailure(Throwable caught) {
					calculation.setCalculated(true);
					Logger.log(caught.getMessage());
				}
			});
		}
	}

	/**
	 * Asks server to start calculation
	 * @param calculation
	 */
	public void startCalculation(final Calculation calculation) {
		Logger.info("Starting calculation for " + calculation.getName());
		calculation.setCalculated(false);
		
		HashMap<String, String> calcData = new HashMap<String, String>();
		for (AmWidget widget : calculation.getInputWidgets()) {
			calcData.put(widget.getDataLink(), widget.getValue());
		}

		calculationService.startCalculation(calcData, config.getWrapperPath()
				+ calculation.getWrapper(),calculation.getName(), new AsyncCallback<Long>() {

			public void onSuccess(Long result) {
				if( result > 0){
					calculation.setTimeID(result);
					Logger.info(calculation.getName() + ": calculation request sent. Id set to: " + result);
				} else {
					Logger.error(calculation.getName() + ": a problem occured, probably wrapper file missing");
					Logger.show();
				}
			}

			public void onFailure(Throwable caught) {
				Logger.error(caught.getMessage());
			}
		});

	}

}
