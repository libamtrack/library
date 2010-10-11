package amtrackservice.client;

import java.util.HashMap;

import com.google.gwt.user.client.rpc.AsyncCallback;

/**
 * The client side stub for the RPC service.
 */
public interface CalculationServiceAsync {
	void startCalculation(HashMap<String, String> input, String wrapperPath, String functionName, AsyncCallback<Long> callback);
	void requestResult(long id, AsyncCallback<HashMap<String,String>> callback);
	void getCalculations(AsyncCallback<HashMap<Long,String>> callback);
}
