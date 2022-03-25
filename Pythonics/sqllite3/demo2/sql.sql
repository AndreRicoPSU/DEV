--CREATE TABLE users (username TEXT NOT NULL PRIMARY KEY ,password NOT NULL, email NOT NULL);

--INSERT INTO users (username, password, email)
--VALUES ("joe", "joepassword", "joe@hahoo.com");

--SELECT * FROM users WHERE username IN ("joe","mike");

--UPDATE users SET email="notmike@gmail.com" WHERE username="mike";

DELETE FROM users WHERE username LIKE "m%e";