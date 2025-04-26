# src/simpyson/gui.py

from PyQt5 import QtWidgets
import sys
from PyQt5.QtWidgets import (
    QFileDialog, QMessageBox, QVBoxLayout, QWidget, QMainWindow, QAction, QInputDialog,
    QSplitter, QListWidget, QHBoxLayout, QDialog, QFormLayout, QLineEdit, QDialogButtonBox
)
from PyQt5.QtCore import Qt
from PyQt5.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go
from simpyson.io import SimpReader
import numpy as np
import os
import copy
from simpyson.converter import hz2ppm, ppm2hz
from simpyson.utils import get_larmor_freq


class SimpysonGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simpyson GUI")
        self.setGeometry(100, 100, 1200, 800)

        self.files_data = {}  # Dictionary to store {filename: data}
        self.current_file = None
        self.data = None
        self.filename = None  # Added to store full path

        self.init_ui()

    def init_ui(self):
        self.create_menu()
        self.create_main_layout()

    def file_list_key_press(self, event):
        if event.key() == Qt.Key_Delete:
            selected_items = self.file_list.selectedItems()
            if selected_items:
                confirmation = QMessageBox.question(
                    self, 'Remove Files', 
                    f'Remove {len(selected_items)} selected file(s)?',
                    QMessageBox.Yes | QMessageBox.No
                )

                if confirmation == QMessageBox.Yes:
                    for item in selected_items:
                        file_name = item.text()

                        if file_name in self.files_data:
                            del self.files_data[file_name]

                        self.file_list.takeItem(self.file_list.row(item))

                    if not self.files_data:
                        self.current_file = None
                        self.data = None
                        self.filename = None
                    elif self.current_file not in self.files_data:
                        self.current_file = next(iter(self.files_data))
                        self.data = self.files_data[self.current_file]['data']
                        self.filename = self.files_data[self.current_file]['path']

                    self.plot_data()
        else:
            QListWidget.keyPressEvent(self.file_list, event)

    def create_main_layout(self):
        # Create central widget with horizontal layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QHBoxLayout(central_widget)

        # Create splitter for resizable panels
        splitter = QSplitter(Qt.Horizontal)

        # Create file list widget with multiple selection
        self.file_list = QListWidget()
        self.file_list.setSelectionMode(QListWidget.ExtendedSelection)
        self.file_list.itemSelectionChanged.connect(self.on_selection_changed)
        splitter.addWidget(self.file_list)
        self.file_list.keyPressEvent = self.file_list_key_press

        # Create plot area
        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        self.browser = QWebEngineView()
        plot_layout.addWidget(self.browser)
        plot_widget.setLayout(plot_layout)

        splitter.addWidget(plot_widget)
        splitter.setSizes([200, 1000])  # Left panel 200px, rest for plot

        layout.addWidget(splitter)

    def create_menu(self):
        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)

        # File Menu
        file_menu = menubar.addMenu('File')

        open_action = QAction('Open', self)
        open_action.triggered.connect(self.open_file)
        open_action.setShortcut('Ctrl+o')
        file_menu.addAction(open_action)

        save_action = QAction('Save', self)
        save_action.triggered.connect(self.save_file)
        save_action.setShortcut('Ctrl+s')
        file_menu.addAction(save_action)

        file_menu.addSeparator()

        exit_action = QAction('Exit', self)
        exit_action.triggered.connect(QtWidgets.qApp.quit)
        file_menu.addAction(exit_action)

        # Process Menu  
        process_menu = menubar.addMenu('Process')

        setup_conversions = QAction('Setup Conversions', self)
        setup_conversions.triggered.connect(self.setup_conversions)
        setup_conversions.setShortcut('Ctrl+g')
        process_menu.addAction(setup_conversions)

        # View Menu
        view_menu = menubar.addMenu('View')

        view_hz = QAction('Hz', self)
        view_hz.triggered.connect(lambda: self.change_view('hz'))
        view_hz.setShortcut('Ctrl+1')
        view_menu.addAction(view_hz)

        view_ppm = QAction('ppm', self)
        view_ppm.triggered.connect(lambda: self.change_view('ppm'))
        view_ppm.setShortcut('Ctrl+2')
        view_menu.addAction(view_ppm)

        view_fid = QAction('FID', self)
        view_fid.triggered.connect(lambda: self.change_view('fid'))
        view_fid.setShortcut('Ctrl+3')
        view_menu.addAction(view_fid)

    def save_file(self):
        if not self.current_file:
            QMessageBox.warning(self, 'Save File', 'No file selected!')
            return

        if not self.data:
            QMessageBox.warning(self, 'Save File', 'No data to save!')
            return

        options = QFileDialog.Options()
        save_filename, _ = QFileDialog.getSaveFileName(
            self,
            'Save File',
            '',
            'All Supported Files (*.csv *.fid *.spe);;CSV Files (*.csv);;FID Files (*.fid);;SPE Files (*.spe)',
            options=options
        )

        if not save_filename: return

        try:
            format = save_filename.lower().split('.')[-1] 

            if format not in ['spe', 'fid', 'csv']:
                QMessageBox.warning(self, 'Save File', 'Unsupported file format!')
                return

            self.data.save(save_filename, format=format)

            QMessageBox.information(self, 'Save File', 'File saved successfully!')
        except Exception as e:
            QMessageBox.warning(self, 'Save File', f'Error saving file: {str(e)}')

    def open_file(self):
        options = QFileDialog.Options()
        filenames, _ = QFileDialog.getOpenFileNames(
            self, 'Open File', '', 'SIMPSON Files (*.spe *.fid)', options=options
        )

        for filename in filenames:
            if filename:
                file_format = filename.split('.')[-1]
                base_name = os.path.basename(filename)
                data = SimpReader(filename, format=file_format)
                view = 'hz' if file_format == 'spe' else 'fid'

                # Store both the data and full path
                self.files_data[base_name] = {
                    'data': data,
                    'path': filename,
                    'view': view,
                }

                # Add to list widget if not already there
                if self.file_list.findItems(base_name, Qt.MatchExactly) == []:
                    self.file_list.addItem(base_name)

                # If this is the first file, display it
                if not self.current_file:
                    self.current_file = base_name
                    self.data = self.files_data[base_name]['data']
                    self.filename = filename
                    self.plot_data()

    def on_selection_changed(self):
        selected_items = self.file_list.selectedItems()
        if selected_items:
            self.current_file = selected_items[0].text()
            self.data = self.files_data[self.current_file]['data']
            self.filename = self.files_data[self.current_file]['path']
            self.plot_data(selected_items)

    def plot_data(self, selected_items=None):
        if not selected_items:
            selected_items = [self.file_list.findItems(self.current_file, Qt.MatchExactly)[0]] if self.current_file else []

        if selected_items:
            fig = go.Figure()

            # Determine x-axis type from first selected item
            first_file = self.files_data[selected_items[0].text()]
            first_data = first_file['data']

            if not 'view' in first_file:
                QMessageBox.warning(self, 'Plot Data', 'Data does not contain valid axis information.')
                return

            match first_file['view']:
                case 'ppm':
                    x_axis = 'ppm'
                    xlabel = 'Chemical Shift (ppm)'
                case 'hz':
                    x_axis = 'hz'
                    xlabel = 'Frequency (Hz)'
                case 'fid':
                    x_axis = 'time'
                    xlabel = 'Time (ms)'

            # Plot each selected spectrum
            for item in selected_items:
                data = self.files_data[item.text()]['data']
                if x_axis in data.data:
                    fig.add_trace(go.Scatter(
                        x=data.data[x_axis],
                        y=data.data['real'],
                        name=item.text(),
                        mode='lines'
                    ))

            # Update layout with consistent axis settings
            fig.update_layout(
                title='NMR Spectrum',
                xaxis_title=xlabel,
                yaxis_title='Intensity',
                showlegend=True
            )

            # Always invert x-axis for ppm and hz
            if x_axis in ['ppm', 'hz']:
                fig.update_xaxes(autorange="reversed")

            # Update plot display
            html_content = fig.to_html(include_plotlyjs='cdn', full_html=True)
            self.browser.setHtml(html_content)
        else:
            self.browser.setHtml('<html><body><h3 style="text-align:center;margin-top:40px;color:#888;">No data to display</h3></body></html>')
            QMessageBox.warning(self, 'Plot Data', 'No data to plot!')

    def change_view(self, view):
        if not (selected_items := self.get_selection()): return
        
        if not all(self.has_setup(item.text()) for item in selected_items):
            QMessageBox.warning(self, 'Not Setup', 'Not all selected items have been setup')
            return
        
        for item in selected_items:
            file_name = item.text()
            file_data = self.files_data[file_name]['data']
            
            view_key = 'time' if view == 'fid' else view

            if not view_key in file_data.data:
                self.convert_to(file_name, view)

            self.files_data[file_name]['view'] = view

        self.plot_data(selected_items)

    def convert_to(self, file_name, target):
        file_data = self.files_data[file_name]['data']

        match target:
            case 'fid':
                if file_data.format == 'fid': return    # Already FID
                new_data = file_data.to_fid()
            case 'hz':
                if file_data.format == 'spe': return    # SPE always has Hz
                new_data = file_data.to_spe()
            case 'ppm':
                if 'ppm' in file_data.data: return      # Already had ppm
                if file_data.format == 'fid':
                    new_data = file_data.to_spe()
                else:
                    new_data = file_data
                hz = new_data.data['hz']
                b0 = new_data.b0
                nucleus = new_data.nucleus
                new_data.data['ppm'] = hz2ppm(hz, b0, nucleus)

        self.files_data[file_name]['data'] = new_data
        if file_name == self.current_file: self.data = new_data

    def setup_conversions(self):
        if not (selected_items := self.get_selection()): return
        
        dialog = QDialog(self)
        dialog.setWindowTitle('Setup Conversions')
        form = QFormLayout(dialog)

        b0_edit = QLineEdit(dialog)
        form.addRow('B0 (e.g, 400MHz or 9.4T:', b0_edit)

        nucleus_edit = QLineEdit(dialog)
        form.addRow('Nucleus (e.g, 1H or 13C:', nucleus_edit)

        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            parent=dialog
        )
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        form.addRow(buttons)

        if dialog.exec_() == QDialog.Accepted:
            b0 = b0_edit.text()
            nucleus = nucleus_edit.text()

            if not (b0 and nucleus):
                return
            
            try:
                _ = get_larmor_freq(b0, nucleus)
            except ValueError as e:
                QMessageBox.warning(self, 'Setup Conversions', f'Invalid Input: {str(e)}')
                return

            for item in selected_items:
                file_name = item.text()
                file_data = self.files_data[file_name]['data']

                file_data.b0 = b0
                file_data.nucleus = nucleus
    
    def get_selection(self):
        selected_items = self.file_list.selectedItems()

        if not selected_items:
            QMessageBox.warning(self, 'No Selection', 'Please select at least one item')
            return None
        
        return selected_items

    def has_setup(self, file_name):
        data = self.files_data[file_name]['data']
        return not (data.b0 is None or data.nucleus is None)

def main():
    app = QtWidgets.QApplication(sys.argv)
    gui = SimpysonGUI()
    gui.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
