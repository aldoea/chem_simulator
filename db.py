import xlrd

sheet = None

def open():
	global sheet
	book=xlrd.open_workbook('db_elements.xlsx')
	sheet=book.sheet_by_index(0)

def getElements():
	global sheet
	components = []
	nrows = sheet.nrows
	for i in range(1,nrows):
		cell = sheet.cell(i,0)
		components.append(cell.value)
	return components	


def getElementValues(row):
    global sheet
    temp = dict()
    a = float(sheet.cell(row, 2).value)
    b = float(sheet.cell(row, 3).value)
    c = float(sheet.cell(row, 4).value)
    tc = float(sheet.cell(row, 5).value)
    pc = float(sheet.cell(row, 6).value)
    temp['A'] = a
    temp['B'] = b
    temp['C'] = c
    temp['Tc'] = tc
    temp['Pc'] = pc
    return temp