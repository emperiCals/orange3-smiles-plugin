import numpy as np
from AnyQt.QtWidgets import QMessageBox, QListWidget, QAbstractItemView
from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.data import Table, Domain, ContinuousVariable
from rdkit import Chem
from rdkit.Chem import Descriptors

class OWSMILESTransformer(OWWidget):
    name = "SMILES Features Transformer"
    description = "Extract molecular descriptors from selected SMILES columns."
    icon = "icons/Unknown.svg"
    priority = 10

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        data = Output("Data", Table)

    want_main_area = False

    def __init__(self):
        super().__init__()
        
        self.data = None

        box = gui.widgetBox(self.controlArea, "SMILES Column Selection")
        
        # 【架构重构】抛弃 Orange 内置的 gui.listBox 绑定机制
        # 强制启用底层原生 Qt 列表组件，以确保数据渲染与传参的绝对可靠
        self.col_listbox = QListWidget()
        self.col_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        box.layout().addWidget(self.col_listbox)

        gui.separator(self.controlArea)
        process_box = gui.widgetBox(self.controlArea, "Execution")
        gui.button(process_box, self, "Process and Output", callback=self.process_and_send)

    @Inputs.data
    def set_data(self, data):
        """接收上游数据并强制重绘视图"""
        self.data = data
        self.col_listbox.clear() # 物理清空现存的所有渲染项

        if data is not None:
            # 遍历所有变量并直接将字符串显式压入 Qt 列表容器中
            for var in data.domain.variables + data.domain.metas:
                self.col_listbox.addItem(var.name)

        self.Outputs.data.send(None)

    def process_and_send(self):
        """执行计算并向下游派发数据"""
        if self.data is None:
            QMessageBox.warning(self, "Error", "No input data available.")
            return

        # 直接从原生 Qt 内存地址中提取用户当前选中的目标节点
        selected_items = self.col_listbox.selectedItems()
        if not selected_items:
            QMessageBox.warning(self, "Error", "Please select at least one SMILES column.")
            return

        # 显式提取文本参数
        selected_col_names = [item.text() for item in selected_items]

        descriptor_names = [
            "MolWt", "ExactMolWt", "HeavyAtomMolWt", "NumValenceElectrons", 
            "FpDensityMorgan1", "FpDensityMorgan2", "FpDensityMorgan3",
            "BertzCT", "Ipc", "HallKierAlpha", "Kappa1", "Kappa2", "Kappa3",
            "Chi0", "Chi1", "Chi0n", "Chi1n", "Chi0v", "Chi1v",
            "MolLogP", "MolMR", "TPSA", "NumHAcceptors", "NumHDonors", "NumRotatableBonds"
        ]

        num_rows = len(self.data)
        num_cols_to_process = len(selected_col_names)
        num_descriptors = len(descriptor_names)

        new_features_matrix = np.zeros((num_rows, num_cols_to_process * num_descriptors))
        new_attributes = list(self.data.domain.attributes)

        col_idx = 0
        for col_name in selected_col_names:
            var = self.data.domain[col_name]

            for desc_name in descriptor_names:
                new_attr_name = f"{col_name}_{desc_name}"
                new_attributes.append(ContinuousVariable(new_attr_name))

            for row_idx in range(num_rows):
                smi_str = str(self.data[row_idx, var]).strip()
                
                if smi_str == '?' or smi_str == 'nan' or smi_str == '':
                    mol = None
                else:
                    mol = Chem.MolFromSmiles(smi_str)

                if mol is not None:
                    new_features_matrix[row_idx, col_idx * num_descriptors + 0] = Descriptors.MolWt(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 1] = Descriptors.ExactMolWt(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 2] = Descriptors.HeavyAtomMolWt(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 3] = Descriptors.NumValenceElectrons(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 4] = Descriptors.FpDensityMorgan1(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 5] = Descriptors.FpDensityMorgan2(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 6] = Descriptors.FpDensityMorgan3(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 7] = Descriptors.BertzCT(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 8] = Descriptors.Ipc(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 9] = Descriptors.HallKierAlpha(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 10] = Descriptors.Kappa1(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 11] = Descriptors.Kappa2(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 12] = Descriptors.Kappa3(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 13] = Descriptors.Chi0(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 14] = Descriptors.Chi1(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 15] = Descriptors.Chi0n(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 16] = Descriptors.Chi1n(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 17] = Descriptors.Chi0v(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 18] = Descriptors.Chi1v(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 19] = Descriptors.MolLogP(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 20] = Descriptors.MolMR(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 21] = Descriptors.TPSA(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 22] = Descriptors.NumHAcceptors(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 23] = Descriptors.NumHDonors(mol)
                    new_features_matrix[row_idx, col_idx * num_descriptors + 24] = Descriptors.NumRotatableBonds(mol)
                else:
                    new_features_matrix[row_idx, col_idx * num_descriptors : (col_idx + 1) * num_descriptors] = np.nan

            col_idx += 1

        new_domain = Domain(new_attributes, self.data.domain.class_vars, self.data.domain.metas)

        if self.data.X is not None and self.data.X.size > 0:
            new_X = np.hstack((self.data.X, new_features_matrix))
        else:
            new_X = new_features_matrix

        out_table = Table.from_numpy(new_domain, X=new_X, Y=self.data.Y, metas=self.data.metas)
        out_table.name = f"{self.data.name}_transformed"

        self.Outputs.data.send(out_table)