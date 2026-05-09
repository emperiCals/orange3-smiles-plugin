import numpy as np
from io import BytesIO
from AnyQt.QtWidgets import QMessageBox, QLabel, QScrollArea, QSizePolicy
from AnyQt.QtGui import QPixmap, QImage
from AnyQt.QtCore import Qt
from Orange.widgets.widget import OWWidget, Input
from Orange.widgets import gui
from Orange.data import Table
from rdkit import Chem
from rdkit.Chem import Draw

class OWSMILESViewer(OWWidget):
    name = "SMILES Viewer"
    description = "Render 2D topological molecular structures from selected SMILES data."
    icon = "icons/Unknown.svg"
    priority = 20

    class Inputs:
        data = Input("Data", Table)

    want_main_area = True

    def __init__(self):
        super().__init__()
        
        self.data = None
        self.available_columns = []
        self.selected_column_index = 0
        self.smiles_list = []
        self.selected_row_index = []

        box_col = gui.widgetBox(self.controlArea, "SMILES Column Selection")
        self.col_combo = gui.comboBox(
            box_col, self, "selected_column_index", items=self.available_columns,
            callback=self.update_smiles_list
        )

        box_row = gui.widgetBox(self.controlArea, "Molecule Instance Selection")
        self.row_listbox = gui.listBox(
            box_row, self, "selected_row_index", "smiles_list",
            selectionMode=gui.listBox.SingleSelection
        )

        gui.separator(self.controlArea)
        process_box = gui.widgetBox(self.controlArea, "Rendering Execution")
        gui.button(process_box, self, "Render Selected Molecule", callback=self.render_molecule)

        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignCenter)
        self.image_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(self.image_label)
        
        self.mainArea.layout().addWidget(scroll_area)

    @Inputs.data
    def set_data(self, data):
        self.data = data
        self.available_columns.clear()
        self.col_combo.clear()
        self.smiles_list.clear()
        self.image_label.clear()
        
        if data is not None:
            for var in data.domain.variables + data.domain.metas:
                self.available_columns.append(var.name)
                self.col_combo.addItem(var.name)
            
            if self.available_columns:
                self.selected_column_index = 0
                self.update_smiles_list()

    def update_smiles_list(self):
        self.smiles_list.clear()
        self.row_listbox.clear()
        if self.data is not None and self.available_columns:
            col_name = self.available_columns[self.selected_column_index]
            var = self.data.domain[col_name]
            for row_idx in range(len(self.data)):
                val = str(self.data[row_idx, var]).strip()
                self.smiles_list.append(f"Row {row_idx + 1}: {val}")
                self.row_listbox.addItem(f"Row {row_idx + 1}: {val}")

    def render_molecule(self):
        if not self.selected_row_index:
            QMessageBox.warning(self, "Error", "No molecular instance selected.")
            return

        row_idx = self.selected_row_index[0]
        col_name = self.available_columns[self.selected_column_index]
        var = self.data.domain[col_name]
        smi_str = str(self.data[row_idx, var]).strip()

        if smi_str == '?' or smi_str == 'nan' or smi_str == '':
            QMessageBox.warning(self, "Error", "Invalid or missing SMILES string.")
            self.image_label.clear()
            return

        mol = Chem.MolFromSmiles(smi_str)
        if mol is None:
            QMessageBox.warning(self, "Error", "RDKit structure resolution failed.")
            self.image_label.clear()
            return

        try:
            img = Draw.MolToImage(mol, size=(600, 600))
            bio = BytesIO()
            img.save(bio, format="PNG")
            qimg = QImage()
            qimg.loadFromData(bio.getvalue())
            
            pixmap = QPixmap.fromImage(qimg)
            self.image_label.setPixmap(pixmap)
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Rendering execution failed: {str(e)}")