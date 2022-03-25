import sqlite3
from employee import Employee

#conn = sqlite3.connect("employee.db")
conn  = sqlite3.connect(':memory:')
c = conn.cursor()

#create a table:
c.execute("""CREATE TABLE employees (first text, last text, pay integer)""")
conn.commit()


def insert_emp(emp):
    with conn: #for conn.commit() after this statement
        c.execute("INSERT INTO employees VALUES (:first, :last, :pay)", {'first': emp.first, 'last': emp.last, 'pay': emp.pay})

def get_emps_by_name(lastname):
    c.execute("SELECT * FROM employees WHERE last=:last", {'last': lastname})
    return c.fetchall()

def update_pay(emp, pay):
    with conn:
        c.execute("""UPDATE employees SET pay = :pay 
                    WHERE first = :first and last = :last """, 
                    {'first': emp.first, 'last': emp.last, 'pay': pay})

def remove_emp(emp):
    with conn:
        c.execute("DELETE from employees WHERE first = :first and last = :last",
                    {'first': emp.first, 'last': emp.last})

#Insert with object
emp_1 = Employee('John', 'Doe', 8000)
emp_2 = Employee('Jane', 'Doe', 9000)

insert_emp(emp_1)
insert_emp(emp_2)

emps = get_emps_by_name('Doe')
print(emps)

update_pay(emp_2, 9500)
remove_emp(emp_1)

emps = get_emps_by_name('Doe')
print(emps)

conn.close()




"""""

# another way:

# c.execute("INSERT INTO employee VALUES ('{}', '{}', {}".format(emp_1.first, emp_1.last, emp_1.pay)) # not good way
c.execute("INSERT INTO employee VALUES (?, ?, ?)", (emp_1.first, emp_1.last, emp_1.pay)) # better way
conn.commit()
c.execute("INSERT INTO employee VALUES (:first, :last, :pay)", {'first': emp_2.first, 'last': emp_2.last, 'pay': emp_2.pay}) # better way
conn.commit()

#Select
c.execute("SELECT * FROM employee WHERE last=?", ('SCHAFER',))
print(c.fetchall())
c.execute("SELECT * FROM employee WHERE last=:last", {'last': 'Doe'})
print(c.fetchall())

#Show datas
# c.fetchone()  # Return first registre
# c.fetchall()  # Return all registers
# c.fetchmany(5) # Return first 5 registers

conn.close()

"""""