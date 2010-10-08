package amtrackservice.server;

public class CalculationData {
	private long id;
	private String filePath;
	private boolean calculated;
	private String functionName;
	
	public CalculationData(long id, String filePath, boolean calculated, String functionName) {
		this.filePath = filePath;
		this.calculated =  calculated;
		this.id = id;
		this.functionName = functionName;
	}

	public long getId() {
		return id;
	}

	public void setId(long id) {
		this.id = id;
	}

	public String getFilePath() {
		return filePath;
	}

	public void setFilePath(String filePath) {
		this.filePath = filePath;
	}

	public boolean isCalculated() {
		return calculated;
	}

	public void setCalculated(boolean calculated) {
		this.calculated = calculated;
	}

	public String getFunctionName() {
		return functionName;
	}

	public void setFunctionName(String functionName) {
		this.functionName = functionName;
	}
	
	
}
