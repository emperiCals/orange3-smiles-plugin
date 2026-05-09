import numpy as np
from io import BytesIO
from AnyQt.QtWidgets import QMessageBox, QLabel, QScrollArea, QSizePolicy, QComboBox, QListWidget, QAbstractItemView
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

        # 【架构重构1】抛弃 gui.comboBox，使用原生 QComboBox 处理列选择
        box_col = gui.widgetBox(self.controlArea, "SMILES Column Selection")
        self.col_combo = QComboBox()
        self.col_combo.currentIndexChanged.connect(self.update_smiles_list)
        box_col.layout().addWidget(self.col_combo)

        # 【架构重构2】抛弃 gui.listBox，使用原生 QListWidget 处理行选择
        box_row = gui.widgetBox(self.controlArea, "Molecule Instance Selection")
        self.row_listbox = QListWidget()
        self.row_listbox.setSelectionMode(QAbstractItemView.SingleSelection)
        self.row_listbox.itemSelectionChanged.connect(self.render_molecule)
        box_row.layout().addWidget(self.row_listbox)

        gui.separator(self.controlArea)
        process_box = gui.widgetBox(self.controlArea, "Rendering Execution")
        gui.button(process_box, self, "Render Selected Molecule", callback=self.render_molecule)

        # 图像渲染主视区
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
        
        # 核心逻辑：在修改列表前暂时屏蔽信号，防止引发意外的越界刷新
        self.col_combo.blockSignals(True)
        self.col_combo.clear()
        self.row_listbox.clear()
        self.image_label.clear()
        
        if data is not None:
            # 移除所有类型判断，无条件接收所有域变量
            for var in data.domain.variables + data.domain.metas:
                self.available_columns.append(var.name)
                self.col_combo.addItem(var.name)
                
        # 恢复信号响应
        self.col_combo.blockSignals(False)
        
        # 数据加载完成后，自动触发一次分子列表的加载
        if self.col_combo.count() > 0:
            self.update_smiles_list()

    def update_smiles_list(self):
        """当用户在下拉菜单切换不同数据列时，刷新下方的分子行列表"""
        self.row_listbox.clear()
        self.image_label.clear()
        
        if self.data is None or self.col_combo.count() == 0:
            return

        col_idx = self.col_combo.currentIndex()
        if col_idx < 0:
            return
            
        col_name = self.available_columns[col_idx]
        var = self.data.domain[col_name]
        
        # 遍历选定列的所有行并显示
        for row_idx in range(len(self.data)):
            val = str(self.data[row_idx, var]).strip()
            self.row_listbox.addItem(f"Row {row_idx + 1}: {val}")
            
        # 默认选中第一行，触发图像渲染
        if self.row_listbox.count() > 0:
            self.row_listbox.setCurrentRow(0)

    def render_molecule(self):
        """当用户选中某一行分子时，执行底层 RDKit 渲染"""
        selected_items = self.row_listbox.selectedItems()
        if not selected_items:
            return

        row_idx = self.row_listbox.currentRow()
        col_idx = self.col_combo.currentIndex()
        
        if row_idx < 0 or col_idx < 0 or self.data is None:
            return
            
        col_name = self.available_columns[col_idx]
        var = self.data.domain[col_name]
        smi_str = str(self.data[row_idx, var]).strip()

        if smi_str == '?' or smi_str == 'nan' or smi_str == '':
            self.image_label.setText("Invalid or missing SMILES string.")
            return

        mol = Chem.MolFromSmiles(smi_str)
        if mol is None:
            self.image_label.setText("Error: RDKit structure resolution failed.")
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
            self.image_label.setText(f"Rendering execution failed: {str(e)}")