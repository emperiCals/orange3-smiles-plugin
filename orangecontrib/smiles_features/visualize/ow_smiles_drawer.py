import numpy as np
from io import BytesIO
from AnyQt.QtWidgets import QMessageBox, QLabel, QScrollArea, QSizePolicy
from AnyQt.QtGui import QPixmap, QImage
from AnyQt.QtCore import Qt
from Orange.widgets.widget import OWWidget
from Orange.widgets import gui
from rdkit import Chem
from rdkit.Chem import Draw

class OWSingleSMILESDrawer(OWWidget):
    name = "Single SMILES Drawer"
    description = "Independently draw a molecular structure by manually inputting a SMILES string."
    icon = "icons/Unknown.svg"
    priority = 30

    # 声明此组件不需要任何前置节点的数据输入
    inputs = []
    outputs = []

    want_main_area = True

    def __init__(self):
        super().__init__()
        
        self.smiles_input_text = ""

        # 构建左侧控制面板区：输入框与执行按钮
        box_input = gui.widgetBox(self.controlArea, "Independent SMILES Input")
        self.smiles_edit = gui.lineEdit(
            box_input, self, "smiles_input_text",
            placeholderText="Type SMILES (e.g., c1ccccc1)",
            callback=self.render_molecule
        )
        
        gui.separator(self.controlArea)
        gui.button(box_input, self, "Render Molecule", callback=self.render_molecule)

        # 构建右侧主视区：图像渲染容器
        self.image_label = QLabel()
        self.image_label.setAlignment(Qt.AlignCenter)
        self.image_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(self.image_label)
        
        self.mainArea.layout().addWidget(scroll_area)

    def render_molecule(self):
        smi_str = self.smiles_input_text.strip()
        
        if not smi_str:
            self.image_label.clear()
            return

        # 调用 RDKit 底层引擎解析 SMILES
        mol = Chem.MolFromSmiles(smi_str)
        if mol is None:
            self.image_label.setText("Error: Invalid SMILES syntax or RDKit parsing failure.")
            return

        # 执行二维拓扑结构绘制
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