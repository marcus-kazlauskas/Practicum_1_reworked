package com.company;

public class TmaException extends Exception{
    /**
     * Throw exception when tridiagonal matrix algorithm doesn't reach desired accuracy
     */
    private final String message;

    public TmaException(String message) {
        this.message = message;
    }

    @Override
    public String getMessage() {
        return message;
    }
}
