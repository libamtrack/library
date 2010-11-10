package amtrackservice.client;

import amtrackservice.client.gui.MainView;

/**
 * This class provides the access to the TextBox on the main GUI, to print out
 * warnings and errors.
 * 
 * @author Christoph Kolb, 2010 Hochschule Heilbronn
 */
public class Logger {

	private static MainView mv = null;

	/**
	 * sets the MainView on which the logs will be printed.
	 * @param mv
	 */
	public static void init(MainView mv) {
		Logger.mv = mv;
	}

	/**
	 * adds a new line on the TextBox of MainView
	 * @param message the message to be logged
	 */
	public static void log(String message) {
		//if (mv != null);
		//	mv.addConsoleLine(message);
	}
	
	/**
	 * adds a new line on the TextBox of MainView
	 * @param message the message to be logged
	 */
	public static void info(String message) {
		//if (mv != null);
		//	mv.addConsoleLine("[INFO]: "+message);
	}
	
	/**
	 * adds a new line on the TextBox of MainView
	 * @param message the message to be logged
	 */
	public static void error(String message) {
		//if (mv != null);
		//	mv.addConsoleLine("[ERROR]: "+message);
	}
}
