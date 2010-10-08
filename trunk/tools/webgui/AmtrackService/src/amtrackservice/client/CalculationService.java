package amtrackservice.client;

import java.util.HashMap;

import com.google.gwt.user.client.rpc.RemoteService;
import com.google.gwt.user.client.rpc.RemoteServiceRelativePath;

/**
 * The client side stub for the RPC service.
 */
@RemoteServiceRelativePath("greet")
public interface CalculationService extends RemoteService {
	long startCalculation(HashMap<String, String> input, String wrapperPath, String functionName) throws IllegalArgumentException;
	HashMap<String,String> requestResult(long id) throws IllegalArgumentException;
	HashMap<Long,String> getCalculations() throws IllegalArgumentException;
}
