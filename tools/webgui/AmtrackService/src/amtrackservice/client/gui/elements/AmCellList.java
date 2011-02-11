package amtrackservice.client.gui.elements;

import java.util.HashMap;
import java.util.List;
import java.util.UUID;

import amtrackservice.client.MapList;

import com.google.gwt.cell.client.ButtonCell;
import com.google.gwt.cell.client.CheckboxCell;
import com.google.gwt.cell.client.EditTextCell;
import com.google.gwt.cell.client.FieldUpdater;
import com.google.gwt.core.client.GWT;
import com.google.gwt.safehtml.shared.SafeHtmlUtils;
import com.google.gwt.uibinder.client.UiBinder;
import com.google.gwt.uibinder.client.UiField;
import com.google.gwt.user.cellview.client.CellTable;
import com.google.gwt.user.cellview.client.Column;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.Widget;
import com.google.gwt.view.client.HasData;
import com.google.gwt.view.client.ListDataProvider;
import com.google.gwt.view.client.MultiSelectionModel;
import com.google.gwt.view.client.SelectionModel;

/**
 * The data source for values information 
 */
class ValuesDatabase {

  /**
   * The provider that holds the list of values in the database.
   */
  private ListDataProvider<Double> dataProvider = new ListDataProvider<Double>();

  
  /**
   * Construct a new contact database.
   */
  ValuesDatabase() {}

  
  /**
   * Add a new value
   * @param value to add.
   */
  public void addValue(Double value) {
    List<Double> values = dataProvider.getList();
    // Remove the contact first so we don't add a duplicate.
    values.remove(value);
    values.add(value);
  }

  /**
   * Remove a value
   * @param value to remove.
   */
  public void removeValue(Double value) {
    List<Double> values = dataProvider.getList();
    values.remove(value);
  }
  
  /**
   * Clears database
   */
  public void clear() {
    List<Double> values = dataProvider.getList();
    values.clear();
  }
  
  /**
   * Get values
   * @return values
   */
  public List<Double> getValues() {
    return dataProvider.getList();
  }
  
  /**
   * Add a display to the database. 
   * @param display a {@Link HasData}.
   */
  public void addDataDisplay(HasData<Double> display) {
    dataProvider.addDataDisplay(display);
  }

  public ListDataProvider<Double> getDataProvider() {
    return dataProvider;
  }

  /**
   * Refresh all displays.
   */
  public void refreshDisplays() {
	  dataProvider.refresh();
  }

}


public class AmCellList extends AmWidget {

	private ValuesDatabase database;
	
	@UiField(provided = true)
	public CellTable<Double> cellTable;
	
	private Widget widget;

	/**
	 * Add the columns to the table.
	 */
	private void initTableColumns( final SelectionModel<Double> selectionModel ) {

		Column<Double, Boolean> checkColumn = new Column<Double, Boolean>( new CheckboxCell(true)) {
			public Boolean getValue(Double object) {
				return selectionModel.isSelected(object);
			}
		};
		cellTable.addColumn(checkColumn, SafeHtmlUtils.fromSafeConstant("<br/>"));

		// Values
		Column<Double, String> valuesColumn = new Column<Double, String>( new EditTextCell()) {
			public String getValue(Double object) {
				return object.toString();
			}
		};

		String valuesColumnName = "Value";
		cellTable.addColumn(valuesColumn, valuesColumnName);
		valuesColumn.setFieldUpdater(new FieldUpdater<Double, String>() {
			public void update(int index, Double object, String value) {
				database.removeValue(object);
				try{
					object = Double.parseDouble(value);
				} catch (NumberFormatException e){
					System.out.println("Problem: " + e.getMessage());
					object = 0.0;
				}
				database.addValue(object);
				database.refreshDisplays();
			}
		});
	
		// Remove buttons
		String removeButtonsColumnName = "Buttons";
		Column<Double, String> removeButtonColumn = new Column<Double, String>(  new ButtonCell()) {
			public String getValue(Double object) {
				return "Remove";
			}
		};
		cellTable.addColumn(removeButtonColumn, removeButtonsColumnName);
		removeButtonColumn.setFieldUpdater(new FieldUpdater<Double, String>() {
			public void update(int index, Double object, String value) {
				database.removeValue(object);
				database.refreshDisplays();
			}
		});

	}
	
	interface Binder extends UiBinder<Widget, AmCellList> {
	}
	
	private String data;
	private MapList<String, String> preset;

	public AmCellList(String label, String datatype, HTML description,
			MapList<String, String> preset, String dataX) {
		super(label, datatype, description);
		
		database = new ValuesDatabase();
		
		data = dataX;
		this.preset = preset;
		
		System.out.println("New AmCellList, data = " + data);
		
		this.cellTable = new CellTable<Double>();
		this.cellTable.setWidth("100%");		

		// Add a selection model so we can select cells.
		final SelectionModel<Double> selectionModel = new MultiSelectionModel<Double>( );
		cellTable.setSelectionModel(selectionModel);

		initTableColumns(selectionModel);

		this.setDefault();
		
		// Add the CellList to the adapter in the database.
		database.addDataDisplay(cellTable);

		Binder uiBinder = GWT.create(Binder.class);
		this.widget = uiBinder.createAndBindUi(this);

	}

	@Override
	public String getValue() {
		String s = "";
		for (Double d : database.getValues()) {
			s += " " + d.toString();
		}
		return s;
	}

	@Override
	public Widget getWidget() {
		return this.widget;
	}

	@Override
	public int setValue(HashMap<String, String> valueMap) {
		database.clear();
		return this.appendValue(valueMap);
	}

	@Override
	public String getDataLink() {
		return data;
	}

	@Override
	public void setDefault() {
		database.clear();
		if( this.preset != null ){
			for (String s : preset.getKeys()) {
				database.addValue(Double.valueOf(preset.getValue(s)));
			}
		}
	}

	@Override
	public int appendValue(HashMap<String, String> valueMap) {
		String str = valueMap.get(data);
		if (str != null) {
			String[] s = str.split(" ");
			for (String c : s) {
				Double d = Double.valueOf(c);
				database.addValue(d);
			}
			System.out.println("Return code = 0, data = " + data + " str = " + str);
			return 0;
		}
		System.out.println("Return code = -1, data = " + data + " str = " + str);
		return -1;
	}

}
