package amtrackservice.client;

import amtrackservice.client.Calculation;

/**
 * This class provides the access to the TextBox on the main GUI, to print out
 * warnings and errors.
 * 
 * @author Christoph Kolb, 2010 Hochschule Heilbronn
 */
public class Logger {

	private static Calculation c = null;

	/**
	 * sets the MainView on which the logs will be printed.
	 * @param mv
	 */
	public static void init(Calculation c) {
		Logger.c = c;
	}

	/**
	 * adds a new line on the TextBox of MainView
	 * @param message the message to be logged
	 */
	public static void log(String message) {
		if( c != null )
			c.addLoggerMessage("[LOG]: " + message);
		System.out.println("[LOG]: " + message);
	}
	
	/**
	 * adds a new line on the TextBox of MainView
	 * @param message the message to be logged
	 */
	public static void info(String message) {
		if( c != null )
			c.addLoggerMessage("[INFO]: " + message);
		System.out.println("[INFO]: " + message);
	}
	
	/**
	 * adds a new line on the TextBox of MainView
	 * @param message the message to be logged
	 */
	public static void error(String message) {
		if( c != null )
			c.addLoggerMessage("[ERROR]: " + message);
		System.out.println("[ERROR]: " + message);
	}
}
